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

#ifndef WIN32 

#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <boost/thread/future.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/array.hpp>

// there is some strange linker error on widows. can't be bothered now... visual c++ will probably do better whole program optimization than gcc anyway...
#define PWDIST_INLINE
#endif
#include "pairwise_seq_distance.h"

#include "ivymike/tdmatrix.h"
#include "ivymike/algorithm.h"
// #include "ivymike/cycle.h"

#include "parsimony.h"
#include "pvec.h"
#include "fasta.h"
#include "ivymike/tree_parser.h"
#include "pars_align_seq.h"

#include <ivymike/time.h>

using namespace std;
using namespace ivy_mike::tree_parser_ms;


template<class pvec_t>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {
//     static int ct;
    //vector<parsimony_state> m_pvec;
  
    vector<uint8_t> m_raw_seq;
    pvec_t m_pvec;
public:
//     int m_ct;
    my_adata_gen() {

//         cout << "my_adata\n";

    }

    virtual ~my_adata_gen() {

//         cout << "~my_adata\n";

    }

    virtual void visit() {
//         cout << "tr: " << m_ct << "\n";
    }
    void init_pvec(const vector< uint8_t >& seq) {
        m_raw_seq = seq;
        reset_pvec();
//        m_pvec.init( seq );
    }
    vector<uint8_t>&get_raw_seq() {
        return m_raw_seq;
    }
    
    void reset_pvec() {
        m_pvec.init( m_raw_seq );
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

// template<class pvec_t>
// class newview_service {
//     boost::thread_group m_tg;
//     boost::barrier m_barrier;
//     boost::barrier m_barrier2;
//     
//     volatile bool m_finish;
//     
//     static void work_outer( size_t rank, newview_service<pvec_t> *this_ ) {
//         this_->work(rank);
//     }
//     
//     void work( size_t rank ) {
//         m_barrier2.wait();
//         while( !m_finish ) {
//             cout << "waiting: " << rank << "\n";
//             m_barrier.wait();
//             
//             
//             cout << "waited: " << rank << "\n";
//             
//             m_barrier2.wait();
//         }
//         cout << "finished: " << rank << "\n";
//     }
//     
// public:
//     newview_service( size_t n_threads ) : m_barrier( n_threads + 1 ), m_barrier2( n_threads + 1 ), m_finish(false) {
//         while( m_tg.size() < n_threads ) {
//             m_tg.create_thread( boost::bind( &newview_service::work_outer, m_tg.size(), this ));
//             
//         }
//         
//         
//   //      m_barrier2.wait();
//         
//     }
//     
//     void do_it() {
//         m_barrier2.wait();
//         m_barrier.wait();
//         
//     }
//     
//     ~newview_service() {
//         m_finish = true;
//         m_barrier2.wait();
//         m_tg.join_all();
//     }
// };

template<class pvec_t>
void do_newview( pvec_t &root_pvec, lnode *n1, lnode *n2, bool incremental ) {
    typedef my_adata_gen<pvec_t> my_adata;

    deque<rooted_bifurcation<lnode> > trav_order;

//     cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";

    rooted_traveral_order( n1, n2, trav_order, incremental );
//     cout << "traversal: " << trav_order.size() << "\n";

    for( deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
//         cout << *it << "\n";

        my_adata *p = dynamic_cast<my_adata *>( it->parent->m_data.get());
        my_adata *c1 = dynamic_cast<my_adata *>( it->child1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( it->child2->m_data.get());
//         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;

//         cout << "tip case: " << (*it) << "\n";
        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->child1->backLen, it->child2->backLen, it->tc);

    }





    {
        my_adata *c1 = dynamic_cast<my_adata *>( n1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( n2->m_data.get());

//         tip_case tc;

        if( c1->isTip && c2->isTip ) {
//                 cout << "root: TIP TIP\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_TIP );
        } else if( c1->isTip && !c2->isTip ) {
//                 cout << "root: TIP INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c2->get_pvec();
        } else if( !c1->isTip && c2->isTip ) {
//                 cout << "root: INNER TIP\n";
            pvec_t::newview(root_pvec, c2->get_pvec(), c1->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c1->get_pvec();
        } else {
//                 cout << "root: INNER INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, INNER_INNER );
        }


    }
//     cout << hex;
//     for( vector< parsimony_state >::const_iterator it = root_pvec.begin(); it != root_pvec.end(); ++it ) {
//         cout << *it;
//     }
//
//     cout << dec << endl;

}


class step_add {
    typedef pvec_cgap pvec_t;
    
    typedef my_adata_gen<pvec_t> my_adata;
    typedef my_fact_gen<my_adata> my_fact;
    //auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    
    const static int score_match = 3;
    const static int score_match_cgap = -4;
    const static int score_gap_open = -3;
    const static int score_gap_extend = -1;
    const static bool align_global = true;
    
    ivy_mike::tree_parser_ms::ln_pool m_ln_pool;
    string m_seq_file_name;
    vector<string> m_qs_names;
    vector<vector<uint8_t> > m_qs_seqs;
//    vector<vector<uint8_t> > m_qs_seqs_mapped; // the content of m_qs_seqs, but mapped by m_pw_scoring_matrix
    
    vector<vector <uint8_t> > m_qs_nongappy;
    
    ivy_mike::tdmatrix<float> m_pw_dist;
    scoring_matrix m_pw_scoring_matrix;
    boost::dynamic_bitset<> m_used_seqs;
    lnode *m_tree_root;
    pars_align_seq::arrays m_seq_arrays;
        
    vector<lnode *> m_leafs;
    align_arrays_traceback<int> m_align_arrays_traceback;
    
    std::ofstream m_inc_ali;
    
    static void seq_to_nongappy_pvec( vector<uint8_t> &seq, vector<uint8_t> &pvec ) {
        pvec.resize( 0 );
        
        for( unsigned int i = 0; i < seq.size(); i++ ) {
            uint8_t ps = dna_parsimony_mapping::d2p(seq[i]);
            
            if( ps == 0x1 || ps == 0x2 || ps == 0x4 || ps == 0x8 ) {
                pvec.push_back(ps);
            }
            
        }
        
    }
    
    
    struct ali_task {
        boost::promise<int> m_prom;
        pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > m_edge;
        const vector<uint8_t> &m_qs_pvec; // WARNINIG: this is most likely a reference to a local object, which is kept in scope until the promise is fullfilled (but not longer)
        
        ali_task( pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > edge, const vector<uint8_t> &qs_pvec ) 
        : m_edge(edge), m_qs_pvec(qs_pvec) {
            
        }
        
        void work() {
            m_prom.set_value(42);
        }
        
    };
    
    
    static void ali_work( step_add *sa, int rank ) {
        pvec_t root_pvec;
        vector<uint8_t> seq_tmp;
        vector<uint8_t> aux_tmp; 
        align_arrays<int32_t> arr;
        
        while( true ) {
            auto_ptr<ali_task> t(sa->get_task());
            
            if( t.get() == 0 ) {
                break;
            }
            
            {
                boost::lock_guard<boost::mutex> tree_lock( sa->m_t_mutex );
                
                do_newview( root_pvec, t->m_edge.first, t->m_edge.second, sa->m_incremental_newview );
                sa->m_incremental_newview = true;
            }
            
            root_pvec.to_int_vec(seq_tmp);
            root_pvec.to_aux_vec(aux_tmp);
            int res2 = align_pvec_score<int32_t,align_global>(seq_tmp, aux_tmp, t->m_qs_pvec, score_match, score_match_cgap, score_gap_open, score_gap_extend, arr );
//             cout << "thread align: " << res2 << "\n";
            
            t->m_prom.set_value(res2);
            
        }
        cout << "worker exit\n";
    }
    
    static void ali_work_vec( step_add *sa, int rank ) {
        // TODO: code review! This went to smoothly, I don' trust it...
        
        const size_t W = 8;
        typedef short score_t;
        
        
        
        
        vector<uint8_t> seq_tmp;
        vector<uint8_t> aux_tmp; 
        
        
        
        vector<ali_task *> tasks;
        vector<pvec_t> root_pvecs(W);
        aligned_buffer<short> a_prof;//(seq_tmp.size() * W);
        aligned_buffer<short> a_aux_prof; //(seq_tmp.size() * W);
        aligned_buffer<short> out_score(W);
        align_vec_arrays<score_t> arr;
        
        while( true ) {
            assert( tasks.empty() );
            
//             auto_ptr<ali_task> t(sa->get_task());
            
            size_t ref_len = 0;
            
//             for( size_t i = 0; i < W; ++i ) {
//                 ali_task *t = sa->get_task().release();
//                 
//                 if( t == 0 ) {
//                     break;
//                 }
//                 
//                 
//                 tasks.push_back(t);
//             }
            
            sa->get_task_block(tasks, W);
            
            if( tasks.empty() ) {
                break;
            }
            
            
//             std::cout << "tasks size: " << tasks.size() << "\n";

            
            
            {
                boost::lock_guard<boost::mutex> tree_lock( sa->m_t_mutex );
            
                for( size_t i = 0; i < tasks.size(); ++i ) {
                    do_newview( root_pvecs[i], tasks[i]->m_edge.first, tasks[i]->m_edge.second, sa->m_incremental_newview );
                    sa->m_incremental_newview = true;
                
                    if( i == 0 ) {
                        ref_len = root_pvecs[i].size();
                    } else {
                        if( ref_len != root_pvecs[i].size() ) {
                            throw std::runtime_error( "quirk: ref_len != root_pvecs[i].size()" );
                        }
                        
                    }
                    
                    if( a_prof.size() < ref_len * W ) {
                        // zero initialize on resize, so that valgrind will not complain about uninitialized reads if tasks.size() < W after resize...
                        a_prof.resize( ref_len * W, 0 );
                        a_aux_prof.resize( ref_len * W, 0 );
                    }
                    
                    // FIXME: to_aux_vec sets a_aux_prof to 0xFFFF for cgap columns. The correct value is dependent on the width of the vector unit.
                    root_pvecs[i].to_int_vec_strided<aligned_buffer<short>::iterator,W>(a_prof.begin() + i);
                    root_pvecs[i].to_aux_vec_strided<aligned_buffer<short>::iterator,W>(a_aux_prof.begin() + i);
                }
                
                
            }
            align_pvec_score_vec<short,W,align_global>(a_prof, a_aux_prof, tasks[0]->m_qs_pvec, score_match, score_match_cgap, score_gap_open, score_gap_extend, out_score, arr );            
            
            sa->m_thread_ncup.at(rank) += uint64_t(ref_len) * tasks.size() * tasks[0]->m_qs_pvec.size();
            
            for( size_t i = 0; i < tasks.size(); ++i ) {
                tasks[i]->m_prom.set_value( out_score[i] );
            }
            
            for( std::vector< ali_task* >::iterator it = tasks.begin(); it != tasks.end(); ++it ) {
                delete *it;
            }
            tasks.clear();

            
        }
        cout << "worker exit\n";
    }
    
    
    auto_ptr<ali_task> get_task() {
        boost::unique_lock<boost::mutex> lock( m_q_mutex );
        
        while( m_queue.empty() && !m_queue_exit ) {
            m_q_cond.wait( lock );
        }
        
        if( m_queue_exit ) {
            return auto_ptr<ali_task>(0);
        } else {
            ali_task *t = m_queue.front();
            m_queue.pop_front();
            return auto_ptr<ali_task>(t);
        }
    }
    
    bool get_task_block( vector<ali_task *> &tasks, size_t max_block_size ) {
        boost::unique_lock<boost::mutex> lock( m_q_mutex );
        
        while( m_queue.empty() && !m_queue_exit ) {
            m_q_cond.wait( lock );
        }
        
        if( m_queue_exit ) {
            return false;
        } else {
            assert( !m_queue.empty() );
            while( tasks.size() < max_block_size && !m_queue.empty() ) {
                ali_task *t = m_queue.front();
                m_queue.pop_front();
                tasks.push_back(t);
//                 std::cout << "task push: " << tasks.size() << "\n";
            }
            
            return true;
        }
    }
    
    // shared ali_task queue: 
    // protected by m_q_mutex
    boost::mutex m_q_mutex;
    deque<ali_task *>m_queue;
    bool m_queue_exit;

    boost::condition_variable m_q_cond;
    // end-of protected by m_q_mutex
    
    // shared tree: 
    // protected by m_t_mutex
    boost::mutex m_t_mutex;
    bool m_incremental_newview;
    // declared above: lnode *m_tree_root;
    // newviews inside the ali_task worker threads must be protected.
    // during modifications in the main thread the workes must be blocked
    
    // end-of protected by m_t_mutex
 
    boost::thread_group m_thread_group;
    
    std::vector<uint64_t> m_thread_ncup;
    
public:
    step_add( const char *seq_name ) 
    : m_ln_pool( new my_fact() ),
    m_seq_file_name(seq_name),
    m_pw_scoring_matrix(3,0),
    m_seq_arrays(true),
    m_queue_exit(false)
    {
        {
            ifstream qsf( m_seq_file_name.c_str() );
            read_fasta( qsf, m_qs_names, m_qs_seqs);
        } 
        m_used_seqs.resize( m_qs_names.size() );
        
        
        while( m_thread_group.size() < 4 ) {
            int rank = int(m_thread_group.size());
            m_thread_group.create_thread( boost::bind( &step_add::ali_work_vec, this, rank ) );
            m_thread_ncup.push_back(0);
        }
        
        if( true ) {
            m_inc_ali.open( "inc_ali.txt" );
        }
        
    }
    
    ~step_add() {
        // cleanup the ali_task worker threads
        m_q_mutex.lock();
        m_queue_exit = true;
        m_q_mutex.unlock();
        m_q_cond.notify_all();
        
        
        m_thread_group.join_all();
        
            
    }


    pair<size_t,size_t> calc_dist_matrix() {
        m_pw_dist.init_size(m_qs_names.size(), m_qs_names.size());
    
        vector<vector<uint8_t> > qs_mapped;
        qs_mapped.reserve(m_qs_names.size() );
        
        // pre-map raw qs seqs to 'state numbers' (=scoring matrix rows/columns)
        for( vector< vector< uint8_t > >::iterator it = m_qs_seqs.begin(); it != m_qs_seqs.end(); ++it) 
        {
            qs_mapped.push_back(vector< uint8_t >());//(it->size()));
            qs_mapped.back().reserve(it->size());
            
            for_each( it->begin(), it->end(), scoring_matrix::valid_state_appender<vector< uint8_t > >(m_pw_scoring_matrix, qs_mapped.back() ));
        }
        
        
        
        ivy_mike::tdmatrix<int> out_scores(m_qs_names.size(), m_qs_names.size());
        
        
        if( false ) {
            ifstream is( "out_scores.bin" );
            
            is.seekg(0, ios_base::end);
            streampos size = is.tellg();
            is.seekg(0, ios_base::beg);
            if( size != out_scores.num_elements() * sizeof(int)) {
                throw runtime_error( "bad external outscores\n" );
            }
            is.read((char*)out_scores.begin(), sizeof(int) * out_scores.num_elements() );
        } else {
            pairwise_seq_distance(qs_mapped, out_scores, m_pw_scoring_matrix, -5, -2, 4);
            ofstream os( "out_scores.bin" );
            os.write((char*)out_scores.begin(), sizeof(int) * out_scores.num_elements() );
        }
        
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
        
        vector<uint8_t> atmp = m_qs_seqs.at(a);
        vector<uint8_t> btmp = m_qs_seqs.at(b);
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
        
        vector<float> dist_sum;
        
        while( f != m_used_seqs.npos ) {
            ivy_mike::odmatrix<float> slice = m_pw_dist[f];
            if( dist_sum.empty() ) {
                dist_sum.assign( slice.begin(), slice.end() );
            } else {
                ivy_mike::binary_twizzle( dist_sum.begin(), dist_sum.end(), slice.begin(), dist_sum.begin(), plus<float>() );
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
    void gapstream_to_alignment( const std::vector<uint8_t> &gaps, const std::vector<uint8_t> &raw, std::vector<uint8_t> &out, uint8_t gap_char, bool upper ) {

        std::vector<uint8_t>::const_reverse_iterator rit = raw.rbegin();


        // 'gap indicator': if upper is set, insert gap into reference (=if gaps[i] == 2).
        const uint8_t gap_ind = upper ? 2 : 1;
        
        
        
        for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {

            
            
            if ( *git == gap_ind ) {
                out.push_back(gap_char);
            } else {
                out.push_back(*rit);
                ++rit;
            } 
        }

        std::reverse(out.begin(), out.end());
    }
    bool insertion_step() {
        size_t candidate = find_next_candidate();
        
        cout << "candidate: " << candidate << "\n";
        
        bool valid = candidate != size_t(-1);
        
        if( !valid ) {
            return false;
        }
        
        
        m_used_seqs[candidate] = true;
        
        
        
//         copy( m_qs_seqs[candidate].begin(), m_qs_seqs[candidate].end(), ostream_iterator<char>(cout) );
//         cout << "\n";
        
        edge_collector<lnode> ec;
        visit_edges( m_tree_root, ec);
        
        cout << "edges: " << ec.m_edges.size() << "\n";
        
        
        
        vector<uint8_t> qs_pvec;
        seq_to_nongappy_pvec(m_qs_seqs.at(candidate), qs_pvec);
        
        //
        // search for best insertion edge
        //
        pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > best_edge;
        int best_score = INT_MIN;
        vector<uint8_t> best_tb;
        
        ivy_mike::timer t1;
        vector<uint8_t> seq_tmp;
        vector<uint8_t> aux_tmp;    
        
        
//         
        const size_t W = 8;
        aligned_buffer<short> out_score(W);
        
        ivy_mike::perf_timer perf_timer;
        
        
        deque<ali_task *> tasks;
        vector<boost::shared_future<int> > futures;
        futures.reserve(ec.m_edges.size());
        for( vector< pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > >::iterator it = ec.m_edges.begin(); it != ec.m_edges.end(); ++it ) {
            tasks.push_back( new ali_task(*it, qs_pvec) );
            
            futures.push_back( boost::shared_future<int>(tasks.back()->m_prom.get_future()) );
        }
        
        {
            boost::lock_guard<boost::mutex> lock( m_q_mutex );
            m_incremental_newview = false;
            m_queue.swap(tasks);
        }
        
        ivy_mike::timer thread_timer;
        m_q_cond.notify_all();
        
        assert( futures.size() == ec.m_edges.size() );
        for( size_t i = 0; i < ec.m_edges.size(); ++i ) {
            int res = futures[i].get();
            //cout << "result: " << i << " = " << res << "\n";
            if( res > best_score ) {

                best_score = res;
                best_edge = ec.m_edges[i];
            }
        }
        //double eticks2 = getticks();
        perf_timer.add_int();
        
        //
        // at this point all worker threads must be blocking on the empty queue (TODO: maybe add explicit check).
        // it is now safe again to motify the tree
        //
        {
            //             cout << "ticks: " << eticks1 << " " << eticks2 << " " << eticks3 << "\n";
            double dt = thread_timer.elapsed();
            size_t sum_ncup = std::accumulate( m_thread_ncup.begin(), m_thread_ncup.end(), uint64_t(0) );
            std::fill( m_thread_ncup.begin(), m_thread_ncup.end(), 0 );
            cout << sum_ncup << " in " << dt << "s : " << sum_ncup / (dt * 1e9) << " GNCUP/s\n";
        }
        
        //
        // generate traceback for best insertion position
        //
        {
            vector<uint8_t> seq_tmp;
            vector<uint8_t> aux_tmp;    
            
            pvec_t root_pvec;
            
            do_newview( root_pvec, best_edge.first, best_edge.second, true );
        
            root_pvec.to_int_vec(seq_tmp);
            root_pvec.to_aux_vec(aux_tmp);
            
            
            
//             size_t stride = 1;
//             size_t aux_stride = 1;
//             pars_align_seq pas( &seq_tmp[0], &qs_pvec[0], seq_tmp.size(), qs_pvec.size(), stride, &aux_tmp[0], aux_stride, m_seq_arrays, 0 );
           // int res = pas.alignFreeshift(INT_MAX);

//             cout << "size: " << seq_tmp.size() << " " << qs_pvec.size() << "\n";
            
            int res2;
            if( align_global ) {
                res2 = align_global_pvec<int>(seq_tmp, aux_tmp, qs_pvec, score_match, score_match_cgap, score_gap_open, score_gap_extend, best_tb, m_align_arrays_traceback );
            } else {
                res2 = align_freeshift_pvec<int>(seq_tmp, aux_tmp, qs_pvec, score_match, score_match_cgap, score_gap_open, score_gap_extend, best_tb, m_align_arrays_traceback );
            }
            
            if( res2 != best_score ) {
                cerr << "res2 != best_score: " << res2 << " " << best_score << "\n";
               throw runtime_error( "alignment error" );
            }
        }
        
        //
        // apply traceback to reference sequences
        //
        
        for( std::vector< ivy_mike::tree_parser_ms::lnode* >::iterator it = m_leafs.begin(); it != m_leafs.end(); ++it ) {
            my_adata *adata = (*it)->m_data->get_as<my_adata>();

            // directly alter the 'raw_seq' stored in the my_adata ojects of the tip-node
            
            vector< uint8_t > &raw_seq = adata->get_raw_seq();
            vector< uint8_t > new_seq;
            // update old raw_seq with new gaps and swap in the new sequence
            gapstream_to_alignment(best_tb, raw_seq, new_seq, '-', true );
            raw_seq.swap(new_seq);
//             cout << "len: " << raw_seq.size() << " " << best_tb.size() << "\n";
            // re-create internal pvec
            adata->reset_pvec();
        }
        
        //
        // apply traceback to query sequence and insert new tip into tree
        //
        
        vector<uint8_t> tip_seq;
        gapstream_to_alignment(best_tb, m_qs_seqs[candidate], tip_seq, '-', false );
//         cout << "tip len: " << tip_seq.size() << "\n";
        
        
        lnode *nc = lnode::create( m_ln_pool );
        nc->m_data->setTipName(m_qs_names.at(candidate));
        nc->m_data->isTip = true;
        
        nc->m_data->get_as<my_adata>()->init_pvec(tip_seq);
        
        lnode *nn = lnode::create( m_ln_pool );
        parser::twiddle(nc, nn, 1.0, "I1", 0 );
        
        
        best_edge.first->back = 0;
        best_edge.second->back = 0;
        parser::twiddle( best_edge.first, nn->next, 1.0, "I1", 0 );
        parser::twiddle( best_edge.second, nn->next->next, 1.0, "I1", 0 );
        
        
        m_leafs.push_back(nc);
        perf_timer.add_int();
#if 0
        {
            stringstream ss;
            ss << "tree_" << setw(5) << setfill('0') << m_leafs.size();
            
            ofstream os( ss.str().c_str() );
//         cout << ">>>>>>>>>>>>>>\n";
            print_newick( next_non_tip( towards_tree( m_tree_root )), os);
//         cout << "<<<<<<<<<<<<<<\n";
        }
        {
            stringstream ss;
            ss << "ali_" << setw(5) << setfill('0') << m_leafs.size();
            
            ofstream os( ss.str().c_str() );
//         cout << ">>>>>>>>>>>>>>\n";
            //print_newick( next_non_tip( towards_tree( m_tree_root )), os);
//         cout << "<<<<<<<<<<<<<<\n";
            
            
            
            for( std::vector< ivy_mike::tree_parser_ms::lnode* >::const_iterator it = m_leafs.begin(); it != m_leafs.end(); ++it ) {
                my_adata *adata = (*it)->m_data->get_as<my_adata>();
                
                copy( adata->get_raw_seq().begin(), adata->get_raw_seq().end(), ostream_iterator<char>(os) );
            
                os << "\n";
            }
            
            
            
        }
#endif   
        
        if( m_inc_ali.good() ) {
            m_inc_ali << m_leafs.size() << "\n";
            for( std::vector< ivy_mike::tree_parser_ms::lnode* >::const_iterator it = m_leafs.begin(); it != m_leafs.end(); ++it ) {
                my_adata *adata = (*it)->m_data->get_as<my_adata>();
                
                copy( adata->get_raw_seq().begin(), adata->get_raw_seq().end(), ostream_iterator<char>(m_inc_ali) );
                
                m_inc_ali << "\n";
            }
            
            m_inc_ali << "\n\n";
            
        }
        
        perf_timer.add_int();
        perf_timer.print();
        
//         std::cout << "ticks: " << eticks2 - eticks1 << " " << eticks3 - eticks2 << " " << eticks4 - eticks3 << "\n";
        
        
        return valid;
           
    }
    void write_phylip( ostream &os ) {
        //
        // write sequences in dfs-ordering, so that topologically close sequences cluster in the output phylip file
        //
        
        tip_collector<lnode>tc;
        
        visit_lnode(m_tree_root, tc );
        
        if( tc.m_nodes.empty() ) {
            throw runtime_error( "tip_collector: empty" );
        }
        
        size_t max_name_len = 0;
        for( vector< tr1::shared_ptr< ivy_mike::tree_parser_ms::lnode > >::const_iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
            max_name_len = max(max_name_len, (*it)->m_data->tipName.size());
        }
        size_t seq_len = tc.m_nodes.front()->m_data->get_as<my_adata>()->get_raw_seq().size();
        os << tc.m_nodes.size() << " " << seq_len << "\n";
        for( vector< tr1::shared_ptr< ivy_mike::tree_parser_ms::lnode > >::const_iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
            my_adata *adata = (*it)->m_data->get_as<my_adata>();

            
            
            os << setw(max_name_len + 1) << left << setfill( ' ' ) << adata->tipName;
            copy( adata->get_raw_seq().begin(), adata->get_raw_seq().end(), ostream_iterator<char>(os) );
            
            os << "\n";
        }
    }
    void write_newick( ostream &os ) {
        print_newick( next_non_tip( towards_tree( m_tree_root )), os); 
    }
    
};

int main( int argc, char **argv ) {
    //mapped_file qsf( "test_1604/1604.fa" );
	
//     string ta = "GATTACAGATTACA";
//     string tb = "GATTACAGATTA";
//     
//     vector<uint8_t> a(ta.begin(), ta.end());
//     vector<uint8_t> b(tb.begin(), tb.end());
//     
//     scoring_matrix sm(3,0);
//     
//     
//   
//     
//     align_freeshift( sm, a, b, -5, -3);
//     return 0;
    
    
//     const char *filename = (argc == 2) ? argv[1] : "test_150/150.fa";
     const char *filename = (argc == 2) ? argv[1] : "test_218/218.fa";

    
    step_add sa(filename);
    pair<size_t,size_t> start_pair = sa.calc_dist_matrix();
    
    cout << "start: " << start_pair.first << " " << start_pair.second << "\n";
    sa.start_tree( start_pair.first, start_pair.second );
    
    
    ivy_mike::timer t1;
    
    
    
    while( sa.insertion_step() ) {
        
    }
    
    {
        ofstream os_ali( "sa_alignment.phy" );
        sa.write_phylip( os_ali );
        
        ofstream os_tree( "sa_tree.phy" );
        sa.write_newick( os_tree );
    }
    
    
    cout << "time: " << t1.elapsed() << "\n";
}
