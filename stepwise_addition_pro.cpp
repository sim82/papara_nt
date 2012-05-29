/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cctype>

#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "pairwise_seq_distance.h"
#include "stepwise_align.h"
#include "raxml_interface.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/smart_ptr.h"
#include "ivymike/tree_parser.h"
#include "ivymike/tdmatrix.h"
#include "ivymike/algorithm.h"
#include "ivymike/tree_traversal_utils.h"


namespace tree_parser = ivy_mike::tree_parser_ms;

using ivy_mike::tree_parser_ms::ln_pool;
using ivy_mike::tree_parser_ms::lnode;

using ivy_mike::apply_lnode;
using ivy_mike::back_insert_ifer;
using ivy_mike::iterate_lnode;
using ivy_mike::rooted_bifurcation;
namespace ublas = boost::numeric::ublas;


typedef std::vector<unsigned char> sequence;

class addition_order {
public:



    addition_order( const scoring_matrix &sm, const std::vector<sequence> &mapped_seqs )
     : used_seqs_( mapped_seqs.size() )
    {
        const size_t num_seqs = mapped_seqs.size();

        std::cout << "size: " << num_seqs << "\n";
        pw_dist_.init_size(num_seqs, num_seqs);


        ivy_mike::tdmatrix<int> out_scores( mapped_seqs.size(), mapped_seqs.size() );

        const size_t num_ali_threads = 4;
        pairwise_seq_distance(mapped_seqs, out_scores, sm, -5, -2, num_ali_threads, 64);
        init_pw_dist_from_msa_score_matrix(out_scores);

    }


    size_t find_next_candidate() {

        if( pw_dist_.size() == 0 ) {
            throw std::runtime_error( "find_next_candidate called with empty pw-dist matrix");
        }

        if( dist_acc_.empty() ) {
            //
            // create initial distance accumulator if it does not exist.
            //
            size_t f = used_seqs_.find_first();

            std::vector<float> dist_sum;

            while( f != used_seqs_.npos ) {
                ivy_mike::odmatrix<float> slice = pw_dist_[f];
                if( dist_sum.empty() ) {
                    dist_sum.assign( slice.begin(), slice.end() );
                } else {
                    std::transform( dist_sum.begin(), dist_sum.end(), slice.begin(), dist_sum.begin(), std::plus<float>() );
                }

                f = used_seqs_.find_next(f);
            }
            dist_acc_.swap( dist_sum );
        }

        float min_dist = 1e8;
        size_t min_element = size_t(-1);
        for( size_t i = 0; i < dist_acc_.size(); i++ ) {
            if( !used_seqs_[i] && dist_acc_[i] < min_dist ) {
                min_dist = dist_acc_[i];
                min_element = i;
            }
        }

        assert( min_element != size_t(-1) || used_seqs_.count() == used_seqs_.size() );

        if( min_element != size_t(-1) ) {

            // update accumulator
            assert( min_element != size_t(-1));
            ivy_mike::odmatrix<float> slice = pw_dist_[min_element];

            assert( slice.size() == dist_acc_.size() );

            // element-wise calculate dist_acc_ = dist_acc_ + slice;
            std::transform( dist_acc_.begin(), dist_acc_.end(), slice.begin(), dist_acc_.begin(), std::plus<float>() );
            used_seqs_[min_element] = true;
        }


        return min_element;

    }


    std::pair<size_t,size_t> first_pair() const {
        return first_pair_;
    }

private:
    void init_pw_dist_from_msa_score_matrix( ivy_mike::tdmatrix<int> &out_scores ) {
        size_t li = -1, lj = -1;
        float lowest_dist = 1e8;
        int min = *(std::min_element( out_scores.begin(), out_scores.end() ));
        int max = *(std::max_element( out_scores.begin(), out_scores.end() ));

        for( size_t i = 0; i < out_scores.size(); i++ ) {

            for( size_t j = 0; j < out_scores[i].size(); j++ ) {

                // three modes for normalizing: min, max and mean
                //const float norm = min( ma[i][i], ma[j][j] );
                //             const float norm = max( ma[i][i], ma[j][j] );
                const float norm = (out_scores[i][j] - min) / float(max-min);


                const float dist = 1.0 - norm;
                pw_dist_[i][j] = dist;

                if( i != j && dist < lowest_dist ) {
                    lowest_dist = dist;
                    li = i;
                    lj = j;
                }

            }

        }

        used_seqs_[li] = used_seqs_[lj] = true;

        first_pair_ = std::make_pair( li, lj );

    }
    ivy_mike::tdmatrix<float> pw_dist_;
    std::vector<float> dist_acc_;
    boost::dynamic_bitset<> used_seqs_;

    std::pair<size_t,size_t> first_pair_;
   // scoring_matrix scoring_matrix_;
};

class sequences {

public:
    sequences( std::istream &is ) : pw_scoring_matrix_( 3, 0 ) {
        assert( is.good() );

        std::cout << "here\n";
        read_fasta( is, names_, seqs_ );

        std::for_each( seqs_.begin(), seqs_.end(), boost::bind( &sequences::normalize_seq, this, _1) );

        std::vector<std::vector<uint8_t> > qs_mapped;
        mapped_seqs_.reserve(names_.size() );

        // pre-map raw qs seqs to 'state numbers' (=scoring matrix rows/columns)
        for( std::vector< std::vector< uint8_t > >::iterator it = seqs_.begin(); it != seqs_.end(); ++it)
        {
            mapped_seqs_.push_back(std::vector< uint8_t >());//(it->size()));
            mapped_seqs_.back().reserve(it->size());

            std::for_each( it->begin(), it->end(), scoring_matrix::valid_state_appender<std::vector< uint8_t > >(pw_scoring_matrix_, mapped_seqs_.back() ));

            std::copy( mapped_seqs_.back().begin(), mapped_seqs_.back().end(), std::ostream_iterator<int>( std::cout, " " ));

            // the raw sequences stored in 'seqs_' filtered for invalid states (e.g., gaps and characters not
            // present in the scoring matrix already. So the members of 'mapped_seqs_' must have the same length.

            assert( mapped_seqs_.back().size() == it->size() );


        }




//        calc_dist_matrix( false );

    }

    const std::vector<sequence> &mapped_seqs() const {
        return mapped_seqs_;
    }
    const std::vector<sequence> &seqs() const {
        return seqs_;
    }


    const sequence &seq_at( size_t i ) const {
        return seqs_.at(i);
    }

    const scoring_matrix &pw_scoring_matrix() const {
        return pw_scoring_matrix_;
    }

    const std::string &name_at( size_t i ) const {
        return names_.at(i);
    }


    size_t clone_seq( size_t i, const std::string &name ) {
        std::vector<uint8_t> seq = seqs_.at(i);
        std::vector<uint8_t> mapped_seq = mapped_seqs_.at(i);

        ivy_mike::push_back_swap(seqs_, seq);
        ivy_mike::push_back_swap(mapped_seqs_, mapped_seq);

//        seqs_.push_back( seq );
//        mapped_seqs_.push_back( mapped_seq );
        names_.push_back(name);

        return seqs_.size() - 1;
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

    std::vector<std::vector<uint8_t> > seqs_;
    std::vector<std::vector<uint8_t> > mapped_seqs_;
    std::vector<std::string> names_;


    scoring_matrix pw_scoring_matrix_;


};


static void make_tip( lnode *n, const std::string &name ) {
    assert( n != 0 );
    assert( n->m_data != 0 );


    // check if the node cn be a valid tip: at least two back pointers must be null.
    size_t null_back = 0;
    if( n->back == 0 ) {
        ++null_back;
    }
    if( n->next->back == 0 ) {
        ++null_back;
    }
    if( n->next->next->back == 0 ) {
        ++null_back;
    }

    assert( null_back >= 2 );

    n->m_data->isTip = true;
    n->m_data->setTipName( name );
}

static bool has_node_label( lnode *n ) {
    assert( n != 0 );
    assert( n->m_data != 0 );
    return !n->m_data->nodeLabel.empty();
}

class tree_builder {
public:

    tree_builder( sequences * const seqs, addition_order * const order, ln_pool * const pool )
     : seqs_(*seqs),
       order_(order),
       pool_(pool)
    {
        // build the initial tree. This is mostly based on black magic.

        std::pair<size_t,size_t> first = order->first_pair();

        size_t seqa = first.first; // confusing, isn't it?
        size_t seqb = first.second;

    //    std::copy( seqs.seq_at( seqa ).begin(), seqs.seq_at( seqa ).end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";
    //    std::copy( seqs.seq_at( seqb ).begin(), seqs.seq_at( seqb ).end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";


        sequence aligned_a = seqs->seq_at(seqa);
        sequence aligned_b = seqs->seq_at(seqb);

        std::string name_clonea = seqs->name_at( seqa ) + "_clone";
        std::string name_cloneb = seqs->name_at( seqb ) + "_clone";

        size_t seqa_clone = seqs->clone_seq( seqa, name_clonea );
        size_t seqb_clone = seqs->clone_seq( seqb, name_cloneb );

        used_seqs_.resize( seqs->seqs().size(), false );
        used_seqs_[seqa] = true;
        used_seqs_[seqa_clone] = true;
        used_seqs_[seqb] = true;
        used_seqs_[seqb_clone] = true;

        aligned_seqs_.resize( seqs->seqs().size() );

    //    std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";
    //    std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(std::cout));
    //    std::cout << "\n";


        lnode *nx = lnode::create( *pool );
        lnode *ny = lnode::create( *pool );
        tree_parser::twiddle_nodes(nx, ny, 1.0, "MOAL", 0 );


        lnode *na1 = lnode::create( *pool );
        lnode *na2 = lnode::create( *pool );

        lnode *nb1 = lnode::create( *pool );
        lnode *nb2 = lnode::create( *pool );

        make_tip( na1, seqs->name_at( seqa ));
        make_tip( na2, name_clonea);
        make_tip( nb1, seqs->name_at( seqb ));
        make_tip( nb2, name_cloneb);


        tree_parser::twiddle_nodes(na1, nx->next, 1.0, "I1", 0 );
        tree_parser::twiddle_nodes(na2, nx->next->next, 1.0, "I2", 0 );
        tree_parser::twiddle_nodes(nb1, ny->next, 1.0, "I3", 0 );
        tree_parser::twiddle_nodes(nb2, ny->next->next, 1.0, "I4", 0 );

        align_freeshift( seqs->pw_scoring_matrix(), aligned_a, aligned_b, -5, -3 );
        assert( aligned_a.size() == aligned_b.size() );

        aligned_seqs_[seqa_clone] = aligned_a;
        aligned_seqs_[seqa].swap( aligned_a );
        aligned_seqs_[seqb_clone] = aligned_b;
        aligned_seqs_[seqb].swap( aligned_b );

        tree_ = nx;


    }


    void print_matrix( const ublas::matrix<double> &m ) {
       // for( ; first != last)
    }


    bool insertion_step() {

        std::vector<ublas::matrix<double> > pvecs;
        write_ali_and_tree_for_raxml();

        lnode *n = generate_marginal_ancestral_state_pvecs( *pool_, "sa_tree", "sa_ali", &pvecs );

        std::cout << "pvecs: " << pvecs.size() << "\n";

        std::cout << n->backLabel << " " << n->next->backLabel << " " << n->next->next->backLabel << "\n";


//        std::deque<rooted_bifurcation<lnode> > to;
//        rooted_traveral_order_rec(n->next, to, false );
//
//        for( std::deque<rooted_bifurcation<lnode> >::iterator it = to.begin(); it != to.end(); ++it ) {
//            std::cout << *it << "\n";
//        }


        std::vector<lnode *> labelled_nodes;
        iterate_lnode(n, back_insert_ifer( labelled_nodes, has_node_label ));

        std::cout << "num labelled: " << labelled_nodes.size() << "\n";


        lnode *tn = labelled_nodes[1];
        std::cout << tn->m_data->nodeLabel << "\n";

        std::deque<rooted_bifurcation<lnode> > to;
        rooted_traveral_order_rec( tn, to, false );

        for( std::deque<rooted_bifurcation<lnode> >::iterator it = to.begin(); it != to.end(); ++it ) {
            std::cout << *it << "\n";
        }



        size_t next_candidate = order_->find_next_candidate();


        return true;
    }

private:

    void write_ali_and_tree_for_raxml() {
        {
            std::ofstream os( "sa_tree" );
            tree_parser::print_newick( tree_, os );
        }



        std::ofstream os( "sa_ali" );


        size_t pos = used_seqs_.find_first();

        os << used_seqs_.count() << " " << aligned_seqs_.at(pos).size() << "\n";

        while( pos != used_seqs_.npos ){
            assert( pos < aligned_seqs_.size() );

            os << seqs_.name_at(pos) << " ";
            std::copy( aligned_seqs_[pos].begin(), aligned_seqs_[pos].end(), std::ostream_iterator<char>(os));
            os << "\n";

            pos = used_seqs_.find_next(pos);
        }

    }


    const sequences &seqs_;
    addition_order * const order_;
    ln_pool * const pool_;
    boost::dynamic_bitset<> used_seqs_;
    lnode * tree_;

    std::vector<sequence> aligned_seqs_;

};

//void insertion_loop( sequences *seqs, addition_order *order, ln_pool * const pool ) {
//
//    std::pair<size_t,size_t> first = order->first_pair();
//
//    size_t seqa = first.first; // confusing, isn't it?
//    size_t seqb = first.second;
//
////    std::copy( seqs.seq_at( seqa ).begin(), seqs.seq_at( seqa ).end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
////    std::copy( seqs.seq_at( seqb ).begin(), seqs.seq_at( seqb ).end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
//
//
//    sequence aligned_a = seqs->seq_at(seqa);
//    sequence aligned_b = seqs->seq_at(seqb);
//
//    std::string name_clonea = seqs->name_at( seqa ) + "_clone";
//    std::string name_cloneb = seqs->name_at( seqb ) + "_clone";
//
//    size_t seqa_clone = seqs->clone_seq( seqa, name_clonea );
//    size_t seqb_clone = seqs->clone_seq( seqb, name_cloneb );
//
//
//
//
////    std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
////    std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(std::cout));
////    std::cout << "\n";
//
//
//    lnode *nx = lnode::create( *pool );
//    lnode *ny = lnode::create( *pool );
//    tree_parser::twiddle_nodes(nx, ny, 1.0, "MOAL", 0 );
//
//
//    lnode *na1 = lnode::create( *pool );
//    lnode *na2 = lnode::create( *pool );
//
//    lnode *nb1 = lnode::create( *pool );
//    lnode *nb2 = lnode::create( *pool );
//
//    make_tip( na1, seqs->name_at( seqa ));
//    make_tip( na2, name_clonea);
//    make_tip( nb1, seqs->name_at( seqb ));
//    make_tip( nb2, name_cloneb);
//
//
//    tree_parser::twiddle_nodes(na1, nx->next, 1.0, "I1", 0 );
//    tree_parser::twiddle_nodes(na2, nx->next->next, 1.0, "I2", 0 );
//    tree_parser::twiddle_nodes(nb1, ny->next, 1.0, "I3", 0 );
//    tree_parser::twiddle_nodes(nb2, ny->next->next, 1.0, "I4", 0 );
//
//    {
//        std::ofstream os( "sa_tree" );
//        tree_parser::print_newick( nx, os );
//    }
//
//    align_freeshift( seqs->pw_scoring_matrix(), aligned_a, aligned_b, -5, -3 );
//    assert( aligned_a.size() == aligned_b.size() );
//    {
//        std::ofstream os( "sa_ali" );
//        os << "4 " << aligned_a.size() << "\n";
//
//
//        os << seqs->name_at( seqa ) << " ";
//        std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//        os << seqs->name_at( seqa_clone ) << " ";
//        std::copy( aligned_a.begin(), aligned_a.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//
//
//
//
//        os << seqs->name_at( seqb ) << " ";
//        std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//        os << seqs->name_at( seqb_clone ) << " ";
//        std::copy( aligned_b.begin(), aligned_b.end(), std::ostream_iterator<char>(os));
//        os << "\n";
//
//
//    }
//
//
////    size_t next;
////    while( (next = order->find_next_candidate()) != size_t(-1)) {
////        std::cout << "next: " << next << "\n";
////        std::copy( seqs->seq_at( next ).begin(), seqs->seq_at( next ).end(), std::ostream_iterator<char>(std::cout));
////        std::cout << "\n";
////
////    }
//}

int main( int argc, char *argv[] ) {
//    {
//        size_t n = 1024 * 1024 * 1024;
//        boost::dynamic_bitset<> bs( n, false );
//
//        for( size_t i = 0; i < n; ++i ) {
//            if( std::rand() < RAND_MAX / 2 ) {
//                bs.set( i, true );
//            }
//        }
//
//
//        ivy_mike::timer t1;
//        size_t x = 0;
//        size_t pos = bs.find_first();
//        while( pos != bs.npos) {
//            x += pos;
//            pos = bs.find_next(pos);
//        }
//
//
//        std::cout << "elapsed: " << t1.elapsed() << "\n";
//
//        return 0;
//    }



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
//    sptr::shared_ptr<ln_pool> pool;//(new ln_pool(std::auto_ptr<node_data_factory>(new my_fact()) ));

    ln_pool pool;



    std::ifstream sis( filename );
    sequences seqs( sis );

    addition_order order( seqs.pw_scoring_matrix(), seqs.mapped_seqs() );

//    insertion_loop( &seqs, &order, &pool );
    tree_builder builder( &seqs, &order, &pool );




    while(true) {
        bool done = builder.insertion_step();
        if( done ) {
            break;
        }
    }

//    builder.write_ali_and_tree_for_raxml();


}
