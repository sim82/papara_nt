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
#include <iomanip>
#ifndef WIN32 
// there is some strange linker error on widows. can't be bothered now... visual c++ will probably do better whole program optimization than gcc anyway...

#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include <boost/thread/future.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/array.hpp>

#define PWDIST_INLINE
#endif
#include "pairwise_seq_distance.h"

#include "ivymike/tdmatrix.h"
#include "ivymike/algorithm.h"
#include "ivymike/cycle.h"

#include "parsimony.h"
#include "pvec.h"
#include "fasta.h"
#include "ivymike/tree_parser.h"
#include "pars_align_seq.h"

using namespace std;
using namespace ivy_mike::tree_parser_ms;
// void pairwise_seq_distance(vector< vector<uint8_t> > &seq_raw, ivy_mike::tdmatrix<int> &, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );
template<typename score_t>
struct align_arrays {
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
};

template<typename score_t>
score_t align_freeshift_pvec_score( vector<uint8_t> &a, vector<uint8_t> &a_aux, const vector<uint8_t> &b, score_t match_score, score_t match_cgap, score_t gap_open, score_t gap_extend, align_arrays<score_t> &arr ) {
    // meeeep: FIXME!
    
    
    
    
    
    if( arr.s.size() < a.size()  ) {
        arr.s.resize( a.size() );
        arr.si.resize( a.size() );
    }
    
    fill( arr.s.begin(), arr.s.end(), 0 );
    fill( arr.si.begin(), arr.si.end(), 0 );
    const score_t SMALL = -32000;

    score_t max_score = SMALL;
    

    

    
    for( size_t ib = 0; ib < b.size(); ib++ ) {
        int bc = b[ib];
        
        score_t last_sl = SMALL;
        score_t last_sc = 0.0;
        score_t last_sdiag = 0.0;
//         cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << endl;
       // sscore_t *qpp_iter = qprofile(m.state_backmap(bc) * W * asize);
        
        score_t * __restrict s_iter = arr.s.base();
        score_t * __restrict si_iter = arr.si.base();
        score_t * __restrict s_end = s_iter + a.size();
        bool lastrow = ib == (b.size() - 1);
        
        //for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {
            
        for( size_t ia = 0; ia < a.size(); ++ia, ++s_iter, ++si_iter ) {  
            //score_t match = sm.get_score( a[ia], bc );
            
            
            score_t ac = a[ia];
            const bool cgap = a_aux[ia] == AUX_CGAP;
            
                // determine match or mis-match according to parsimony bits coming from the tree.
            score_t match = ( ac & bc ) != 0 ? match_score : 0;
            

            score_t sm = last_sdiag + match;
            
            last_sdiag = *s_iter;

           // score_t sm_zero = sm;

            
            score_t last_sc_OPEN; 
            score_t sl_score_stay;
            
            if( cgap ) {
                last_sc_OPEN = last_sc; 
                sl_score_stay = last_sl;
                sm += match_cgap;
            } else {
                last_sc_OPEN = last_sc + gap_open; 
                sl_score_stay = last_sl + gap_extend;
            }
            
            score_t sl;
            if( sl_score_stay > last_sc_OPEN ) {
                sl = sl_score_stay;
            } else {
                sl = last_sc_OPEN;
            }
            //score_t sl = max( last_sl + gap_extend, last_sc_OPEN );
            
            
            
//             vu::println( GAP_EXT_vec, pb.m_ptr ); getchar();
            
            
            last_sl = sl;
            

            score_t su_gap_open = last_sdiag + gap_open;
            score_t su_GAP_EXTEND = *si_iter + gap_extend;
            
            
            
            score_t su;// = max( su_GAP_EXTEND,  );
            if( su_GAP_EXTEND > su_gap_open ) {
                su = su_GAP_EXTEND;
            } else {
                su = su_gap_open;
            }
            
            
            *si_iter = su;
            
            
            
            
            score_t sc = std::max( sm, std::max( su, sl ) );
            
            last_sc = sc;
            *s_iter = sc;
            
            
            if( s_iter == (s_end - 1) || lastrow ) {
                if( sc  > max_score ) {
                    
                    max_score = sc;
                }
                
                
            }
        }
        
    }
    
    

    return max_score;
    
}

template<typename score_t, size_t W>
void align_freeshift_pvec_score_vec( aligned_buffer<score_t> &a_prof, aligned_buffer<score_t> &a_aux_prof, vector<uint8_t> &b, const score_t match_score, const score_t match_cgap, const score_t gap_open, const score_t gap_extend, aligned_buffer<score_t> &out ) {
    
    
    typedef vector_unit<score_t,W> vu;
    typedef typename vu::vec_t vec_t;
    
    
    // meeeep: FIXME!
    static aligned_buffer<score_t> s;
    static aligned_buffer<score_t> si;
    
    size_t av_size = a_prof.size();
    
    if( s.size() < av_size  ) {
        s.resize( av_size );
        si.resize( av_size );
    }
    
    fill( s.begin(), s.end(), 0 );
    fill( si.begin(), si.end(), 0 );
    const score_t SMALL = vu::SMALL_VALUE;

    vec_t max_score = vu::set1(SMALL);
//     const vec_t full_mask = vu::set1(0xFFFF);
    
    for( size_t ib = 0; ib < b.size(); ib++ ) {
        //int bc = b[ib];
        vec_t bc = vu::set1(b[ib]);
        
        vec_t last_sl = vu::set1(SMALL);
        vec_t last_sc = vu::set1(0);
        vec_t last_sdiag = vu::set1(0);
//         cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << endl;
       // sscore_t *qpp_iter = qprofile(m.state_backmap(bc) * W * asize);
        
        
        score_t * __restrict a_prof_iter = a_prof.base();
        score_t * __restrict a_aux_prof_iter = a_aux_prof.base();
        score_t * __restrict s_iter = s.base();
        score_t * __restrict si_iter = si.base();
        score_t * __restrict s_end = s_iter + av_size;
        bool lastrow = ib == (b.size() - 1);
        
//         vec_t lastrow_mask = vu::set1( lastrow ? 0xFFFF : 0 );
        
        //for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {
            
        //for( size_t ia = 0; ia < a.size(); ++ia, s_iter+=W, si_iter+=W ) {  
        vec_t row_max_score = vu::set1(SMALL);
        bool break_aloop = false;
        for(; !break_aloop; a_prof_iter += W, a_aux_prof_iter += W, s_iter += W, si_iter += W ) {
            //score_t match = sm.get_score( a[ia], bc );
            //vec_t lastcol_mask = (s_iter == (s_end - W)) ? full_mask : lastrow_mask;
            
            break_aloop = s_iter == (s_end - W);
            //vec_t lastcol_mask = vu::bit_or( lastrow_mask, vu::cmp_eq(vu::set1(ia), vu::setzero() ));
            
            
            const vec_t ac = vu::load( a_prof_iter );
            const vec_t cgap = vu::load( a_aux_prof_iter );
            
            
            
            //score_t match = ( ac & bc ) != 0 ? match_score : 0;
            const vec_t non_match = vu::cmp_eq( vu::bit_and( ac, bc ), vu::setzero() );
            
            vec_t match_score_vec = vu::bit_andnot( non_match, vu::set1(match_score));
            const vec_t match_cgap_score = vu::bit_and( cgap, vu::set1(match_cgap));
            match_score_vec = vu::add( match_score_vec, match_cgap_score );
            
            const vec_t gap_open_score = vu::bit_andnot( cgap, vu::set1(gap_open));
            const vec_t gap_extend_score = vu::bit_andnot( cgap, vu::set1(gap_extend));
            
            const vec_t sm = vu::add(last_sdiag, match_score_vec );
            
            //last_sdiag = *s_iter;
            last_sdiag = vu::load( s_iter );
 
            
            const vec_t last_sc_OPEN = vu::add( last_sc, gap_open_score ); 
            const vec_t sl_score_stay = vu::add( last_sl, gap_extend_score );
            
//             if( cgap ) {
//                 last_sc_OPEN = last_sc; 
//                 sl_score_stay = last_sl;
//                 sm += match_cgap;
//             } else {
//                 last_sc_OPEN = last_sc + gap_open; 
//                 sl_score_stay = last_sl + gap_extend;
//             }
            
            const vec_t sl = vu::max( last_sc_OPEN, sl_score_stay );
//             if( sl_score_stay > last_sc_OPEN ) {
//                 sl = sl_score_stay;
//             } else {
//                 sl = last_sc_OPEN;
//             }
            //score_t sl = max( last_sl + gap_extend, last_sc_OPEN );
            
            
            
//             vu::println( GAP_EXT_vec, pb.m_ptr ); getchar();
            
            
            last_sl = sl;
            

            const vec_t su_gap_open = vu::add( last_sdiag, vu::set1(gap_open));
            
            const vec_t si = vu::load( si_iter );
            const vec_t su_GAP_EXTEND = vu::add( si, vu::set1(gap_extend) );
            
            
            
            const vec_t su = vu::max( su_GAP_EXTEND, su_gap_open );// = max( su_GAP_EXTEND,  );
//             if( su_GAP_EXTEND > su_gap_open ) {
//                 su = su_GAP_EXTEND;
//             } else {
//                 su = su_gap_open;
//             }
            
            vu::store( su, si_iter );
            //*si_iter = su;
            
            
            
            
            const vec_t sc = vu::max( sm, vu::max( su, sl ) );
            
            last_sc = sc;
            
            vu::store( sc, s_iter );
//             *s_iter = sc;
            
            
//             vec_t mask = (break_aloop || lastrow) ? full_mask : vu::setzero();
//             max_score = vu::max( max_score, vu::bit_and( sc, mask));
//             if( break_aloop || lastrow ) {
// //                 if( sc  > max_score ) {
// //                     
// //                     max_score = sc;
// //                 }
//                 max_score = vu::max( max_score, sc );
//                 
//             }
            
            row_max_score = vu::max( row_max_score, sc );
        }
	max_score = vu::max( max_score, last_sc );
        // CONTINUE HERE!
        if( lastrow ) {
            max_score = vu::max( max_score, row_max_score );	
	} 
        
        
    }
    
    
    vu::store( max_score, out.base() );
//     return max_score;
    
}


template<typename score_t>
score_t align_freeshift_pvec( vector<uint8_t> &a, vector<uint8_t> &a_aux, vector<uint8_t> &b, score_t match_score, score_t match_cgap, score_t gap_open, score_t gap_extend, vector<uint8_t>& tb_out ) {
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
    
    
    
    if( s.size() < a.size()  ) {
        s.resize( a.size() );
        si.resize( a.size() );
    }
    
    fill( s.begin(), s.end(), 0 );
    fill( si.begin(), si.end(), 0 );
    const score_t SMALL = -32000;

    score_t max_score = SMALL;
    size_t max_a = 0;
    size_t max_b = 0;
    
    vector<bool> sl_stay( a.size() * b.size() );
    vector<bool> su_stay( a.size() * b.size() );
    vector<bool> s_l( a.size() * b.size() );
    vector<bool> s_u( a.size() * b.size() );
    
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
        int bc = b[ib];
        
        score_t last_sl = SMALL;
        score_t last_sc = 0.0;
        score_t last_sdiag = 0.0;
//         cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << endl;
       // sscore_t *qpp_iter = qprofile(m.state_backmap(bc) * W * asize);
        
        score_t * __restrict s_iter = s.base();
        score_t * __restrict si_iter = si.base();
        score_t * __restrict s_end = s_iter + a.size();
        bool lastrow = ib == (b.size() - 1);
        
        //for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {
            
        for( size_t ia = 0; ia < a.size(); ++ia, ++s_iter, ++si_iter ) {  
            //score_t match = sm.get_score( a[ia], bc );
            
            int ac = a[ia];
            const bool cgap = a_aux[ia]  == AUX_CGAP;
            
                // determine match or mis-match according to parsimony bits coming from the tree.
            score_t match = ( ac & bc ) != 0 ? match_score : 0;
            

            score_t sm = last_sdiag + match;
            
            last_sdiag = *s_iter;

            

            
            score_t last_sc_OPEN; 
            score_t sl_score_stay;
            
            if( cgap ) {
                last_sc_OPEN = last_sc; 
                sl_score_stay = last_sl;
                sm += match_cgap;
            } else {
                last_sc_OPEN = last_sc + gap_open; 
                sl_score_stay = last_sl + gap_extend;
            }
            
            score_t sl;
            if( sl_score_stay > last_sc_OPEN ) {
                sl = sl_score_stay;
                sl_stay[ic(ia,ib)] = true;
            } else {
                sl = last_sc_OPEN;
            }
            //score_t sl = max( last_sl + gap_extend, last_sc_OPEN );
            
            
            
//             vu::println( GAP_EXT_vec, pb.m_ptr ); getchar();
            
            
            last_sl = sl;
            

            score_t su_gap_open = last_sdiag + gap_open;
            score_t su_GAP_EXTEND = *si_iter + gap_extend;
            
            
            
            score_t su;// = max( su_GAP_EXTEND,  );
            if( su_GAP_EXTEND > su_gap_open ) {
                su = su_GAP_EXTEND;
                su_stay[ic(ia,ib)] = true;
            } else {
                su = su_gap_open;
            }
            
            
            *si_iter = su;
            
            score_t sc;
            if( (su > sl) && su > sm ) {
                sc = su;
                s_u[ic(ia,ib)] = true;
            } else if( ( sl >= su ) && sl > sm ) {
                sc = sl;
                s_l[ic(ia,ib)] = true;
            } else { // implicit: sm_zero > sl && sm_zero > su
                sc = sm;
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
    
    cout << "max score: " << max_score << "\n";
    cout << "max " << max_a << " " << max_b << "\n";
    
    
    
    
    int ia = a.size() - 1;
    int ib = b.size() - 1;
    
    
    
    assert( ia == max_a || ib == max_b );
    
    bool in_l = false;
    bool in_u = false;
    
    while( ia > max_a ) {
        //         a_tb.push_back(a[ia]);
        //         b_tb.push_back('-');
        tb_out.push_back(1);

        --ia;
//         in_u = true;
    }
    
    while( ib > max_b ) {
        //         a_tb.push_back('-');
        //         b_tb.push_back(b[ib]);
        tb_out.push_back(2);
        --ib;
//         in_l = true;
    }
    
//     cout << "in " << in_l << " " << in_u << "\n";
    
    while( ia >= 0 && ib >= 0 ) {
        size_t c = ic( ia, ib );
        
        if( !in_l && !in_u ) {
            in_l = s_l[c];
            in_u = s_u[c];
            
            if( !in_l && !in_u ) {
                //                 a_tb.push_back(a[ia]);
                //                 b_tb.push_back(b[ib]);
                tb_out.push_back(0);
                --ia;
                --ib;
            }
            
        }
        
        if( in_u ) {
            //             a_tb.push_back('-');
            //             b_tb.push_back(b[ib]);
            tb_out.push_back(2);
            --ib;
            
            in_u = su_stay[c];
        } else if( in_l ) {
            //             a_tb.push_back(a[ia]);
            //             b_tb.push_back('-');
            tb_out.push_back(1);
            --ia;
            
            in_l = sl_stay[c];
        }
        
        
    }
    
    while( ia >= 0 ) {
        //         a_tb.push_back(a[ia]);
        //         b_tb.push_back('-');
        tb_out.push_back(1);
        --ia;
    }
    
    while( ib >= 0 ) {
        //         a_tb.push_back('-');
        //         b_tb.push_back(b[ib]);
        tb_out.push_back(2);
        --ib;
    }
    //reverse( tb_out.begin(), tb_out.end() );
    
//     getchar();
    return max_score;
}


void align_freeshift( const scoring_matrix &sm, vector<uint8_t> &a, vector<uint8_t> &b, float gap_open, float gap_extend ) {

 
    typedef float score_t;
    
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
    
    
    
    if( s.size() < a.size()  ) {
        s.resize( a.size() );
        si.resize( a.size() );
    }
    
    fill( s.begin(), s.end(), 0 );
    fill( si.begin(), si.end(), 0 );
    const score_t SMALL = -1e8;

    score_t max_score = SMALL;
    size_t max_a = 0;
    size_t max_b = 0;
    
    vector<bool> sl_stay( a.size() * b.size() );
    vector<bool> su_stay( a.size() * b.size() );
    vector<bool> s_l( a.size() * b.size() );
    vector<bool> s_u( a.size() * b.size() );
    
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
//         cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << endl;
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
            //score_t sl = max( last_sl + gap_extend, last_sc_OPEN );
            
            
            
//             vu::println( GAP_EXT_vec, pb.m_ptr ); getchar();
            
            
            last_sl = sl;
            

            score_t su_gap_open = last_sdiag + gap_open;
            score_t su_GAP_EXTEND = *si_iter + gap_extend;
            
            
            
            score_t su;// = max( su_GAP_EXTEND,  );
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
    
    cout << "max score: " << max_score << "\n";
    cout << "max " << max_a << " " << max_b << "\n";
    
    
    vector<uint8_t> a_tb;
    vector<uint8_t> b_tb;
    
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
    
//     cout << "in " << in_l << " " << in_u << "\n";
    
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
    copy( a_tb.rbegin(), a_tb.rend(), a.begin() );
    
    b.resize( b_tb.size() );
    copy( b_tb.rbegin(), b_tb.rend(), b.begin() );
    
    copy( a.begin(), a.end(), ostream_iterator<char>(cout) );
    cout << "\n";
    copy( b.begin(), b.end(), ostream_iterator<char>(cout) );
    cout << "\n";
}


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
    
    
    static void ali_work( step_add *sa ) {
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
            int res2 = align_freeshift_pvec_score<int32_t>(seq_tmp, aux_tmp, t->m_qs_pvec, 3, -10, -3, -1, arr );
//             cout << "thread align: " << res2 << "\n";
            
            t->m_prom.set_value(res2);
            
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
    
    deque<ali_task *>m_queue;
    bool m_queue_exit;
    
    boost::mutex m_q_mutex;
    boost::condition_variable m_q_cond;
    boost::mutex m_t_mutex;
    bool m_incremental_newview;
    class queue_cleanup {
        step_add &sa;
    public:
        queue_cleanup( step_add &sa_ ) : sa(sa_) {}
        ~queue_cleanup() {
            sa.m_q_mutex.lock();
            sa.m_queue_exit = true;
            sa.m_q_mutex.unlock();
            sa.m_q_cond.notify_all();
        }
    };
    
    queue_cleanup m_queue_cleanup;
    boost::thread_group m_thread_group;
    
public:
    step_add( const char *seq_name ) 
    : m_ln_pool( new my_fact() ),
    m_seq_file_name(seq_name),
    m_pw_scoring_matrix(3,0),
    m_seq_arrays(true),
    m_queue_exit(false),
    m_queue_cleanup(*this)
    
    {
        {
            ifstream qsf( m_seq_file_name.c_str() );
            read_fasta( qsf, m_qs_names, m_qs_seqs);
        } 
        m_used_seqs.resize( m_qs_names.size() );
        
        
        while( m_thread_group.size() < 4 ) {
            m_thread_group.create_thread( boost::bind( &step_add::ali_work, this ) );
        }
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
        
        
        bool first_edge = true;
        vector<uint8_t> qs_pvec;
        seq_to_nongappy_pvec(m_qs_seqs.at(candidate), qs_pvec);
        
        //
        // search for best insertion edge
        //
        pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > best_edge;
        int best_score = INT_MIN;
        vector<uint8_t> best_tb;
        uint64_t ncup = 0;
        ivy_mike::timer t1;
        vector<uint8_t> seq_tmp;
        vector<uint8_t> aux_tmp;    
        
        double eticks1 = 0.0;
        double eticks2 = 0.0;
        double eticks3 = 0.0;
        
        const size_t W = 8;
        aligned_buffer<short> out_score(W);
        
        
        deque<ali_task *> tasks;
        vector<boost::shared_future<int> > futures;
        futures.reserve(ec.m_edges.size());
        for( vector< pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > >::iterator it = ec.m_edges.begin(); it != ec.m_edges.end(); ++it ) {
//             vector<uint8_t> best_tb2;
            //    cout << it->first->m_data->m_serial << " " << it->second->m_data->m_serial << "\n";
            
            
//             boost::promise<int> prom;
            tasks.push_back( new ali_task(*it, qs_pvec) );
            
            futures.push_back( boost::shared_future<int>(tasks.back()->m_prom.get_future()) );
            
            
#if 0            
            ticks ticks1 = getticks();
            pvec_t root_pvec;
           // first_edge = false;
            do_newview( root_pvec, it->first, it->second, !first_edge );
        
            root_pvec.to_int_vec(seq_tmp);
            root_pvec.to_aux_vec(aux_tmp);

            ticks ticks2 = getticks();
            //int res2 = align_freeshift_pvec_score<int32_t>(seq_tmp, aux_tmp, qs_pvec, 3, -10, -3, -1 );
            ticks ticks3 = getticks();
#if 1
            aligned_buffer<short> a_prof(seq_tmp.size() * W);
            aligned_buffer<short> a_aux_prof(seq_tmp.size() * W);
           for( size_t i = 0; i < seq_tmp.size(); i++ ) {
                fill( a_prof(i * W), a_prof(i*W) + W, seq_tmp[i] );
                short cgap = aux_tmp[i] == AUX_CGAP ? 0xFFFF : 0x0;
                fill( a_aux_prof(i * W), a_aux_prof(i*W) + W, cgap );
                
            }
            align_freeshift_pvec_score_vec<short,W>(a_prof, a_aux_prof, qs_pvec, 3, -10, -3, -1, out_score );
            int res2 = out_score[0];
            for( size_t i = 0; i < W; i++ ) {
                if( res2 != out_score[i] ) {
                    cout << "bad score: " << res2 << " " << out_score[i] << "\n";
                    throw runtime_error( "vector error" );
                }
            }


#endif
            ticks ticks4 = getticks();
            if( res2 > best_score ) {

                best_score = res2;
                best_edge = *it;
                //best_tb = best_tb2;
                //            }
                //
            }
            
            
            ncup += seq_tmp.size() * qs_pvec.size();
            
            eticks1 += elapsed(ticks2, ticks1);
            eticks2 += elapsed(ticks3, ticks2);
            eticks3 += elapsed(ticks4, ticks3);
            first_edge = false;
#endif            
        }
        
        {
            boost::lock_guard<boost::mutex> lock( m_q_mutex );
            m_incremental_newview = false;
            m_queue.swap(tasks);
        }
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
        
        
        {
            cout << "ticks: " << eticks1 << " " << eticks2 << " " << eticks3 << "\n";
            
            cout << ncup << " in " << t1.elapsed() << "s : " << ncup / (t1.elapsed() * 1e9) << " GNCUP/s\n";
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
            
            int res2 = align_freeshift_pvec<short>(seq_tmp, aux_tmp, qs_pvec, 3, -10, -3, -1, best_tb );
            
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
};

int main( int argc, char **argv ) {
    //mapped_file qsf( "test_1604/1604.fa" );
	
    string ta = "GATTACAGATTACA";
    string tb = "GATTACAGATTA";
    
    vector<uint8_t> a(ta.begin(), ta.end());
    vector<uint8_t> b(tb.begin(), tb.end());
    
    scoring_matrix sm(3,0);
    
    
  
    
    align_freeshift( sm, a, b, -5, -3);
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
        ofstream os( "sa_result.phy" );
        sa.write_phylip( os );
    }
    
    cout << "time: " << t1.elapsed() << "\n";
}
