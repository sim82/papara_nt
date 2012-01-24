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
#include <cctype>

#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "tree_utils.h"

#include "stepwise_align.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/smart_ptr.h"
#include "ivymike/tree_parser.h"
#include "ivymike/tdmatrix.h"

using ivy_mike::tree_parser_ms::ln_pool;


class sequences {

public:
    sequences( std::istream &is ) : pw_scoring_matrix_( 3, 0 ) {
        assert( is.good() );

        std::cout << "here\n";
        read_fasta( is, names_, seqs_ );

        std::for_each( seqs_.begin(), seqs_.end(), boost::bind( &sequences::normalize_seq, this, _1) );

        calc_dist_matrix( false );

    }

private:

    void normalize_seq( std::vector<uint8_t> &seq ) {
        std::vector<uint8_t> nseq;
      //  nseq.reserve(seq.size());

        for( std::vector<uint8_t>::iterator it = seq.begin(); it != seq.end(); ++it ) {
            uint8_t c = std::toupper(*it);

            if( pw_scoring_matrix_.state_valid(c)) {
                nseq.push_back(c);
            }
        }

        // shrink-to-fit into original vector 'seq'
       // std::vector<uint8_t> tmp( nseq.begin(), nseq.end() );
        seq.swap(nseq);

    }

    void calc_dist_matrix( bool load_scores ) {


        std::cout << "size: " << names_.size() << "\n";
        pw_dist_.init_size(names_.size(), names_.size());

        std::vector<std::vector<uint8_t> > qs_mapped;
        mapped_seqs_.reserve(names_.size() );

        // pre-map raw qs seqs to 'state numbers' (=scoring matrix rows/columns)
        for( std::vector< std::vector< uint8_t > >::iterator it = seqs_.begin(); it != seqs_.end(); ++it)
        {
            mapped_seqs_.push_back(std::vector< uint8_t >());//(it->size()));
            mapped_seqs_.back().reserve(it->size());

            //std::for_each( it->begin(), it->end(), scoring_matrix::valid_state_appender<std::vector< uint8_t > >(m_pw_scoring_matrix, qs_mapped.back() ));


            std::copy( it->begin(), it->end(),
                    back_insert_ifer(mapped_seqs_.back(), boost::bind(&scoring_matrix::state_valid, &pw_scoring_matrix_, _1)));


       //     std::cout << "mapped: " << mapped_seqs_.back().size() << " " << it->size() << "\n";
        }



//        ivy_mike::tdmatrix<int> out_scores(m_qs_names.size(), m_qs_names.size());
//
//
//        if( load_scores ) {
//            ifstream is( "out_scores.bin" );
//
//            is.seekg(0, ios_base::end);
//            size_t size = is.tellg();
//            is.seekg(0, ios_base::beg);
//            if( size != out_scores.num_elements() * sizeof(int)) {
//                std::cerr << size << " vs " << out_scores.num_elements() * sizeof(int) << "\n";
//
//                throw std::runtime_error( "bad external outscores\n" );
//            }
//            is.read((char*)out_scores.begin(), sizeof(int) * out_scores.num_elements() );
//        } else {
//            pairwise_seq_distance(qs_mapped, out_scores, m_pw_scoring_matrix, -5, -2, m_num_ali_threads, 64);
//            ofstream os( "out_scores.bin" );
//            os.write((char*)out_scores.begin(), sizeof(int) * out_scores.num_elements() );
//        }
//
//        return init_pw_dist_from_local_score_matrix(out_scores);
    }

    std::vector<std::vector<uint8_t> > seqs_;
    std::vector<std::vector<uint8_t> > mapped_seqs_;
    std::vector<std::string> names_;

    ivy_mike::tdmatrix<float> pw_dist_;


    scoring_matrix pw_scoring_matrix_;
};

int main( int argc, char *argv[] ) {
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



    std::map<std::string, std::vector<uint8_t> >out_msa1;
    sptr::shared_ptr<ln_pool> pool;//(new ln_pool(std::auto_ptr<node_data_factory>(new my_fact()) ));



    std::ifstream sis( filename );
    sequences seqs( sis );



}
