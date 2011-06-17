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

#include <functional>

#include <boost/dynamic_bitset.hpp>

#define PWDIST_INLINE
#include "pairwise_seq_distance.h"

#include "ivymike/tdmatrix.h"
#include "ivymike/algorithm.h"

#include "parsimony.h"
#include "pvec.h"
#include "fasta.h"
#include "ivymike/tree_parser.h"
#include "pars_align_seq.h"

using namespace ivy_mike::tree_parser_ms;
// void pairwise_seq_distance(std::vector< std::vector<uint8_t> > &seq_raw, ivy_mike::tdmatrix<int> &, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );


void align_freeshift( const scoring_matrix &sm, std::vector<uint8_t> &a, std::vector<uint8_t> &b, float gap_open, float gap_extend ) {

 
    typedef float score_t;
    
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
    
    
    
    if( s.size() < a.size()  ) {
        s.resize( a.size() );
        si.resize( a.size() );
    }
    
    std::fill( s.begin(), s.end(), 0 );
    std::fill( si.begin(), si.end(), 0 );
    const score_t SMALL = -1e8;

    score_t max_score = SMALL;
    size_t max_a = 0;
    size_t max_b = 0;
    
    std::vector<bool> sl_stay( a.size() * b.size() );
    std::vector<bool> su_stay( a.size() * b.size() );
    std::vector<bool> s_l( a.size() * b.size() );
    std::vector<bool> s_u( a.size() * b.size() );
    
    struct index_calc {
        const size_t as, bs;
        index_calc( size_t as_, size_t bs_ ) : as(as_), bs(bs_) {}
        
        size_t operator()(size_t ia, size_t ib ) {
            assert( ia < as );
            assert( ib < bs );
            
            //return ia * bs + ib;
            return ib * as + ia;
        }
        
    };
    
    index_calc ic( a.size(), b.size() );
    
    
    for( size_t ib = 0; ib < b.size(); ib++ ) {
        char bc = b[ib];
        
        score_t last_sl = SMALL;
        score_t last_sc = 0.0;
        score_t last_sdiag = 0.0;
//         std::cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << std::endl;
       // sscore_t *qpp_iter = qprofile(m.state_backmap(bc) * W * asize);
        
        score_t * __restrict s_iter = s.base();
        score_t * __restrict si_iter = si.base();
        score_t * __restrict s_end = s_iter + a.size();
        bool lastrow = ib == (b.size() - 1);
        
        //for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {
        for( size_t ia = 0; ia < a.size(); ++ia, ++s_iter, ++si_iter ) {  
            score_t match = sm.get_score( a[ia], bc );
            
//             getchar();
            score_t sm = last_sdiag + match;
            
            last_sdiag = *s_iter;

            score_t sm_zero = sm;

            
            score_t last_sc_OPEN = last_sc + gap_open; 
            
            
            score_t sl_score_stay = last_sl + gap_extend;
            score_t sl;
            if( sl_score_stay > last_sc_OPEN ) {
                sl = sl_score_stay;
                sl_stay[ic(ia,ib)] = true;
            } else {
                sl = last_sc_OPEN;
            }
            //score_t sl = std::max( last_sl + gap_extend, last_sc_OPEN );
            
            
            
//             vu::println( GAP_EXT_vec, pb.m_ptr ); getchar();
            
            
            last_sl = sl;
            

            score_t su_gap_open = last_sdiag + gap_open;
            score_t su_GAP_EXTEND = *si_iter + gap_extend;
            
            
            
            score_t su;// = std::max( su_GAP_EXTEND,  );
            if( su_GAP_EXTEND > su_gap_open ) {
                su = su_GAP_EXTEND;
                su_stay[ic(ia,ib)] = true;
            } else {
                su = su_gap_open;
            }
            
            
            *si_iter = su;
            
            score_t sc;
            if( (su > sl) && su > sm_zero ) {
                sc = su;
                s_u[ic(ia,ib)] = true;
            } else if( ( sl >= su ) && sl > sm_zero ) {
                sc = sl;
                s_l[ic(ia,ib)] = true;
            } else { // implicit: sm_zero > sl && sm_zero > su
                sc = sm_zero;
            }

            
            last_sc = sc;
            *s_iter = sc;
            
            
            if( s_iter == s_end - 1 || lastrow ) {
                if( sc > max_score ) {
                    max_a = ia;
                    max_b = ib;
                    max_score = sc;
                }
                
                
            }
        }
        
    }
    
    std::cout << "max score: " << max_score << "\n";
    std::cout << "max " << max_a << " " << max_b << "\n";
    
    
    std::vector<uint8_t> a_tb;
    std::vector<uint8_t> b_tb;
    
    int ia = a.size() - 1;
    int ib = b.size() - 1;
    
    
    
    assert( ia == max_a || ib == max_b );
    
    bool in_l = false;
    bool in_u = false;
    
    while( ia > max_a ) {
        a_tb.push_back(a[ia]);
        b_tb.push_back('-');
        --ia;
//         in_u = true;
    }
    
    while( ib > max_b ) {
        a_tb.push_back('-');
        b_tb.push_back(b[ib]);
        --ib;
//         in_l = true;
    }
    
    std::cout << "in " << in_l << " " << in_u << "\n";
    
    while( ia >= 0 && ib >= 0 ) {
        size_t c = ic( ia, ib );
        
        if( !in_l && !in_u ) {
            in_l = s_l[c];
            in_u = s_u[c];
            
            if( !in_l && !in_u ) {
                a_tb.push_back(a[ia]);
                b_tb.push_back(b[ib]);
                --ia;
                --ib;
            }
            
        }
        
        if( in_u ) {
            a_tb.push_back('-');
            b_tb.push_back(b[ib]);
            --ib;
            
            in_u = su_stay[c];
        } else if( in_l ) {
            a_tb.push_back(a[ia]);
            b_tb.push_back('-');
            --ia;
            
            in_l = sl_stay[c];
        }
        
        
    }
    
     while( ia >= 0 ) {
        a_tb.push_back(a[ia]);
        b_tb.push_back('-');
        --ia;
    }
    
    while( ib >= 0 ) {
        a_tb.push_back('-');
        b_tb.push_back(b[ib]);
        --ib;
    }
    
    a.resize( a_tb.size() );
    std::copy( a_tb.rbegin(), a_tb.rend(), a.begin() );
    
    b.resize( b_tb.size() );
    std::copy( b_tb.rbegin(), b_tb.rend(), b.begin() );
    
    std::copy( a.begin(), a.end(), std::ostream_iterator<char>(std::cout) );
    std::cout << "\n";
    std::copy( b.begin(), b.end(), std::ostream_iterator<char>(std::cout) );
    std::cout << "\n";
}


template<class pvec_t>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {
//     static int ct;
    //std::vector<parsimony_state> m_pvec;
    pvec_t m_pvec;

public:
//     int m_ct;
    my_adata_gen() {

//         std::cout << "my_adata\n";

    }

    virtual ~my_adata_gen() {

//         std::cout << "~my_adata\n";

    }

    virtual void visit() {
//         std::cout << "tr: " << m_ct << "\n";
    }
    void init_pvec(const std::vector< uint8_t >& seq) {


        m_pvec.init( seq );

    }
    pvec_t &get_pvec() {
        return m_pvec;
    }


};

template<class ndata_t>
class my_fact_gen : public ivy_mike::tree_parser_ms::node_data_factory {

    virtual ndata_t *alloc_adata() {

        return new ndata_t;
    }

};

template<class pvec_t>
void do_newview( pvec_t &root_pvec, lnode *n1, lnode *n2, bool incremental ) {
    typedef my_adata_gen<pvec_t> my_adata;

    std::deque<rooted_bifurcation<lnode> > trav_order;

    std::cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";

    rooted_traveral_order( n1, n2, trav_order, incremental );
//     std::cout << "traversal: " << trav_order.size() << "\n";

    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
//         std::cout << *it << "\n";

        my_adata *p = dynamic_cast<my_adata *>( it->parent->m_data.get());
        my_adata *c1 = dynamic_cast<my_adata *>( it->child1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( it->child2->m_data.get());
//         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;

//         std::cout << "tip case: " << (*it) << "\n";
        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->child1->backLen, it->child2->backLen, it->tc);

    }





    {
        my_adata *c1 = dynamic_cast<my_adata *>( n1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( n2->m_data.get());

//         tip_case tc;

        if( c1->isTip && c2->isTip ) {
//                 std::cout << "root: TIP TIP\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_TIP );
        } else if( c1->isTip && !c2->isTip ) {
//                 std::cout << "root: TIP INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c2->get_pvec();
        } else if( !c1->isTip && c2->isTip ) {
//                 std::cout << "root: INNER TIP\n";
            pvec_t::newview(root_pvec, c2->get_pvec(), c1->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c1->get_pvec();
        } else {
//                 std::cout << "root: INNER INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, INNER_INNER );
        }


    }
//     std::cout << std::hex;
//     for( std::vector< parsimony_state >::const_iterator it = root_pvec.begin(); it != root_pvec.end(); ++it ) {
//         std::cout << *it;
//     }
//
//     std::cout << std::dec << std::endl;

}


class step_add {
    typedef pvec_cgap pvec_t;
    
    typedef my_adata_gen<pvec_t> my_adata;
    typedef my_fact_gen<my_adata> my_fact;
    //std::auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    ivy_mike::tree_parser_ms::ln_pool m_ln_pool;
    std::string m_seq_file_name;
    std::vector<std::string> m_qs_names;
    std::vector<std::vector<uint8_t> > m_qs_seqs;
//    std::vector<std::vector<uint8_t> > m_qs_seqs_mapped; // the content of m_qs_seqs, but mapped by m_pw_scoring_matrix
    
    std::vector<std::vector <uint8_t> > m_qs_nongappy;
    
    ivy_mike::tdmatrix<float> m_pw_dist;
    scoring_matrix m_pw_scoring_matrix;
    boost::dynamic_bitset<> m_used_seqs;
    lnode *m_tree_root;
    pars_align_seq::arrays m_seq_arrays;
    
    std::vector<lnode *> m_leafs;
    static void seq_to_nongappy_pvec( std::vector<uint8_t> &seq, std::vector<uint8_t> &pvec ) {
        pvec.resize( 0 );
        
        for( unsigned int i = 0; i < seq.size(); i++ ) {
            uint8_t ps = dna_parsimony_mapping::d2p(seq[i]);
            
            if( ps == 0x1 || ps == 0x2 || ps == 0x4 || ps == 0x8 ) {
                pvec.push_back(ps);
            }
            
        }
        
    }
    
public:
    step_add( const char *seq_name ) 
    : m_ln_pool( new my_fact() ),
    m_seq_file_name(seq_name),
    m_pw_scoring_matrix(3,0),
    m_seq_arrays(true)
    
    {
        {
            std::ifstream qsf( m_seq_file_name.c_str() );
            read_fasta( qsf, m_qs_names, m_qs_seqs);
        } 
        m_used_seqs.resize( m_qs_names.size() );
        
    }

    std::pair<size_t,size_t> calc_dist_matrix() {
        m_pw_dist.init_size(m_qs_names.size(), m_qs_names.size());
    
        std::vector<std::vector<uint8_t> > qs_mapped;
        qs_mapped.reserve(m_qs_names.size() );
        
        // pre-map raw qs seqs to 'state numbers' (=scoring matrix rows/columns)
        for( std::vector< std::vector< uint8_t > >::iterator it = m_qs_seqs.begin(); it != m_qs_seqs.end(); ++it) 
        {
            qs_mapped.push_back(std::vector< uint8_t >());//(it->size()));
            qs_mapped.back().reserve(it->size());
            
            std::for_each( it->begin(), it->end(), scoring_matrix::valid_state_appender<std::vector< uint8_t > >(m_pw_scoring_matrix, qs_mapped.back() ));
        }
        
        
        
        ivy_mike::tdmatrix<int> out_scores(m_qs_names.size(), m_qs_names.size());
        
        pairwise_seq_distance(qs_mapped, out_scores, m_pw_scoring_matrix, -5, -2, 4);
     
        size_t li, lj;
        float lowest_dist = 1e8;
        
        
        for( size_t i = 0; i < out_scores.size(); i++ ) {
            
            for( size_t j = 0; j < out_scores[i].size(); j++ ) {
            
                // three modes for normalizing: min, max and mean
                //const float norm = std::min( ma[i][i], ma[j][j] );
                //             const float norm = std::max( ma[i][i], ma[j][j] );
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
        
        return std::pair<size_t,size_t>(li,lj);
    }
    void start_tree( size_t a, size_t b ) {
        lnode *na = lnode::create( m_ln_pool );
        lnode *nb = lnode::create( m_ln_pool );
        
        m_leafs.push_back(na);
        m_leafs.push_back(nb);
        
        m_used_seqs[a] = true;
        m_used_seqs[b] = true;
        
        assert( ivy_mike::isa<my_adata>( na->m_data.get() ) );
        assert( ivy_mike::isa<my_adata>( nb->m_data.get() ) );
        
        my_adata *da = na->m_data.get()->get_as<my_adata>();
        my_adata *db = nb->m_data.get()->get_as<my_adata>();
        
        std::vector<uint8_t> atmp = m_qs_seqs.at(a);
        std::vector<uint8_t> btmp = m_qs_seqs.at(b);
        scoring_matrix sm( 3, 0 );
        
        align_freeshift(sm, atmp, btmp, -5, -3 );
        
        da->init_pvec(atmp);
        db->init_pvec(btmp);
        
        parser::twiddle(na, nb, 1.0, "I1", 0 );
        na->m_data->setTipName(m_qs_names.at(a));
        nb->m_data->setTipName(m_qs_names.at(b));
        na->m_data->isTip = true;
        nb->m_data->isTip = true;
        
        
        
        
        m_tree_root = na;
    }
    
    size_t find_next_candidate() {
        size_t f = m_used_seqs.find_first();
        
        std::vector<float> dist_sum;
        
        while( f != m_used_seqs.npos ) {
            ivy_mike::odmatrix<float> slice = m_pw_dist[f];
            if( dist_sum.empty() ) {
                dist_sum.assign( slice.begin(), slice.end() );
            } else {
                ivy_mike::binary_twizzle( dist_sum.begin(), dist_sum.end(), slice.begin(), dist_sum.begin(), std::plus<float>() );
            }
            
            f = m_used_seqs.find_next(f);
        }
        
        float min_dist = 1e8;
        size_t min_element = size_t(-1);
        for( size_t i = 0; i < dist_sum.size(); i++ ) {
            if( !m_used_seqs[i] && dist_sum[i] < min_dist ) {
                min_dist = dist_sum[i];
                min_element = i;
            }
        }
        
        return min_element;
    }
    
    bool insertion_step() {
        size_t candidate = find_next_candidate();
        
        std::cout << "candidate: " << candidate << "\n";
        
        bool valid = candidate != size_t(-1);
        
        if( valid ) {
            m_used_seqs[candidate] = true;
        }
        
        
        edge_collector<lnode> ec;
        visit_edges( m_tree_root, ec);
        
        std::cout << "edges: " << ec.m_edges.size() << "\n";
        
        
        bool first_edge = true;
        std::vector<uint8_t> qs_pvec;
        seq_to_nongappy_pvec(m_qs_seqs.at(candidate), qs_pvec);
        for( std::vector< std::pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > >::iterator it = ec.m_edges.begin(); it != ec.m_edges.end(); ++it ) {
            //    std::cout << it->first->m_data->m_serial << " " << it->second->m_data->m_serial << "\n";
            std::vector<int> seq_tmp;
            std::vector<unsigned int> aux_tmp;    
            
            pvec_t root_pvec;
            
            do_newview( root_pvec, it->first, it->second, !first_edge );
        
            root_pvec.to_int_vec(seq_tmp);
            root_pvec.to_aux_vec(aux_tmp);
            
            
            
            size_t stride = 1;
            size_t aux_stride = 1;
            pars_align_seq pas( &seq_tmp[0], &qs_pvec[0], seq_tmp.size(), qs_pvec.size(), stride, &aux_tmp[0], aux_stride, m_seq_arrays, 0 );
            int res = pas.alignFreeshift(INT_MAX);
            
            first_edge = false;
        }
        
        
        return valid;
           
    }
    
};

int main( int argc, char **argv ) {
    //mapped_file qsf( "test_1604/1604.fa" );
	
    std::string ta = "GATTACAGATTACA";
    std::string tb = "GATTACAGATTA";
    
    std::vector<uint8_t> a(ta.begin(), ta.end());
    std::vector<uint8_t> b(tb.begin(), tb.end());
    
    scoring_matrix sm(3,0);
    
    
    align_freeshift( sm, a, b, -5, -3);
//     return 0;
    
    
    const char *filename = (argc == 2) ? argv[1] : "test_1604/1604.fa.100";
    

    
    step_add sa(filename);
    std::pair<size_t,size_t> start_pair = sa.calc_dist_matrix();
    
    std::cout << "start: " << start_pair.first << " " << start_pair.second << "\n";
    sa.start_tree( start_pair.first, start_pair.second );
    
    
    ivy_mike::timer t1;
    
    while( sa.insertion_step() ) {
        
    }
    
    std::cout << "time: " << t1.elapsed() << "\n";
}