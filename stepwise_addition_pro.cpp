/*
 * Copyright (C) 2011 Simon A. Berger
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */
#include "stepwise_align.h"
#include <functional>
#include <iomanip>
#include <numeric>
#include <deque>


#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
//#include <boost/thread/future.hpp>
#include <boost/thread/barrier.hpp>

#include <boost/array.hpp>

#ifndef WIN32
// there is some strange linker error on widows. can't be bothered now... visual c++ will probably do better whole program optimization than gcc anyway...
#define PWDIST_INLINE
#endif
#include "pairwise_seq_distance.h"

#include "ivymike/tdmatrix.h"
#include "ivymike/algorithm.h"
// #include "ivymike/cycle.h"
#include "tree_utils.h"
#include "parsimony.h"
#include "pvec.h"
#include "fasta.h"
#include "ivymike/tree_parser.h"
#include "pars_align_seq.h"
#include "tree_similarity.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/multiple_alignment.h"

#include <boost/tr1/unordered_set.hpp>
#include "ivymike/concurrent.h"
#include "raxml_interface.h"


using std::string;
using std::vector;
using std::map;
using std::pair;
using std::deque;

using std::ios_base;
using std::ifstream;
using std::ofstream;


using namespace ivy_mike::tree_parser_ms;


class queries {


public:

    queries( std::istream &is )
     : m_pw_scoring_matrix(3,0)
    {
        read_fasta( is, names_, seqs_ );

        start_pair_ = calc_dist_matrix(false);
    }


    const ivy_mike::tdmatrix<float> pw_dist() const {
        return m_pw_dist;
    }

    const std::pair<size_t, size_t> start_pair() const {
        return start_pair_;
    }

private:

    std::pair<size_t,size_t> calc_dist_matrix( bool load_scores ) {
        m_pw_dist.init_size(names_.size(), names_.size());

        vector<vector<uint8_t> > qs_mapped;
        qs_mapped.reserve(names_.size() );

        // pre-map raw qs seqs to 'state numbers' (=scoring matrix rows/columns)
        for( vector< vector< uint8_t > >::iterator it = seqs_.begin(); it != seqs_.end(); ++it)
        {
            qs_mapped.push_back(vector< uint8_t >());//(it->size()));
            qs_mapped.back().reserve(it->size());

            for_each( it->begin(), it->end(), scoring_matrix::valid_state_appender<vector< uint8_t > >(m_pw_scoring_matrix, qs_mapped.back() ));
        }



        ivy_mike::tdmatrix<int> out_scores(names_.size(), names_.size());


        if( load_scores ) {
            ifstream is( "out_scores.bin" );

            is.seekg(0, ios_base::end);
            size_t size = is.tellg();
            is.seekg(0, ios_base::beg);
            if( size != out_scores.num_elements() * sizeof(int)) {
                std::cerr << size << " vs " << out_scores.num_elements() * sizeof(int) << "\n";

                throw std::runtime_error( "bad external outscores\n" );
            }
            is.read((char*)out_scores.begin(), sizeof(int) * out_scores.num_elements() );
        } else {
            const size_t num_threads = 2;

            pairwise_seq_distance(qs_mapped, out_scores, m_pw_scoring_matrix, -5, -2, num_threads, 64);
            ofstream os( "out_scores.bin" );
            os.write((char*)out_scores.begin(), sizeof(int) * out_scores.num_elements() );
        }

        return init_pw_dist_from_local_score_matrix(out_scores);
    }

    pair<size_t,size_t> init_pw_dist_from_local_score_matrix( ivy_mike::tdmatrix<int> &out_scores ) {
        size_t li = -1, lj = -1;
        float lowest_dist = 1e8;


        for( size_t i = 0; i < out_scores.size(); i++ ) {

            for( size_t j = 0; j < out_scores[i].size(); j++ ) {

                // three modes for normalizing: min, max and mean
                //const float norm = min( ma[i][i], ma[j][j] );
                //             const float norm = max( ma[i][i], ma[j][j] );
                const float norm = (out_scores[i][i] + out_scores[j][j]) * 0.5;


                int mae;
                if( i <= j ) {
                    mae = out_scores[i][j];
                    //                 mae = ma[j][i];
                } else {
                    mae = out_scores[j][i];

                }

                const float dist = 1.0 - (mae / norm);
                m_pw_dist[i][j] = dist;

                if( i != j && dist < lowest_dist ) {
                    lowest_dist = dist;
                    li = i;
                    lj = j;
                }

            }

        }

        return pair<size_t,size_t>(li,lj);
    }


    std::vector<std::string> names_;
    std::vector<std::vector<uint8_t> > seqs_;

    ivy_mike::tdmatrix<float> m_pw_dist;
    std::pair<size_t,size_t> start_pair_;
    scoring_matrix m_pw_scoring_matrix;



};


class inc_tree {
public:


    inc_tree( queries &qs )
     : qs_(qs)
    {


    }
    size_t find_next_candidate() {

        if( qs_.pw_dist().size() == 0 ) {
            throw std::runtime_error( "find_next_candidate called with empty pw-dist matrix");
        }

        if( m_dist_acc.empty() ) {
            //
            // create initial distance accumulator if it does not exist.
            //
            size_t f = m_used_seqs.find_first();

            vector<float> dist_sum;

            while( f != m_used_seqs.npos ) {
                ivy_mike::odmatrix<float> slice = qs_.pw_dist()[f];
                if( dist_sum.empty() ) {
                    dist_sum.assign( slice.begin(), slice.end() );
                } else {
                    ivy_mike::binary_twizzle( dist_sum.begin(), dist_sum.end(), slice.begin(), dist_sum.begin(), std::plus<float>() );
                }

                f = m_used_seqs.find_next(f);
            }
            m_dist_acc.swap( dist_sum );
        }

        float min_dist = 1e8;
        size_t min_element = size_t(-1);
        for( size_t i = 0; i < m_dist_acc.size(); i++ ) {
            if( !m_used_seqs[i] && m_dist_acc[i] < min_dist ) {
                min_dist = m_dist_acc[i];
                min_element = i;
            }
        }

        assert( min_element != size_t(-1) || m_used_seqs.count() == m_used_seqs.size() );

        if( min_element != size_t(-1) ) {

            // update accumulator
            assert( min_element != size_t(-1));
            ivy_mike::odmatrix<float> slice = qs_.pw_dist()[min_element];

            assert( slice.size() == m_dist_acc.size() );
            ivy_mike::binary_twizzle( m_dist_acc.begin(), m_dist_acc.end(), slice.begin(), m_dist_acc.begin(), std::plus<float>() );
        }
        return min_element;

    }

private:
    const queries &qs_;

    std::vector<float> m_dist_acc;
    boost::dynamic_bitset<> m_used_seqs;
};

int main( int argc, char **argv ) {
	ivy_mike::getopt::parser igp;
	std::string opt_seq_file;

	int num_cores = boost::thread::hardware_concurrency();

	int opt_num_ali_threads;
	int opt_num_nv_threads;
	bool opt_load_scores;

	igp.add_opt('h', false );
	igp.add_opt('f', ivy_mike::getopt::value<std::string>(opt_seq_file) );
	igp.add_opt('j', ivy_mike::getopt::value<int>(opt_num_ali_threads).set_default(num_cores) );
	igp.add_opt('k', ivy_mike::getopt::value<int>(opt_num_nv_threads).set_default(1) );
	igp.add_opt('l', ivy_mike::getopt::value<bool>(opt_load_scores, true).set_default(false) );
	bool ret = igp.parse(argc, argv);

	if( igp.opt_count('h') != 0 || !ret ) {
		std::cout <<
		"  -h        print help message\n";
		return 0;

	}


	if( igp.opt_count('f') != 1 ) {
		std::cerr << "missing option -f\n";
#ifndef WIN32 // hack. make it easier to start inside visual studio
		return 0;
#endif
		opt_seq_file = "test_218/218.fa";
	}

    const char *filename = opt_seq_file.c_str();


    std::ifstream is( opt_seq_file );
    queries qs( is );


}
