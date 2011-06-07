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

#ifndef __align_vec_h
#define __align_vec_h

#include "vec_unit.h"
#include "aligned_buffer.h"
#include "fasta.h"

template <class score_t>
struct persistent_state {
    aligned_buffer<score_t> out;
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;  
};

//#define LIKELY(x) __builtin_expect((x),1)
// #define LIKELY(x) (x)
template <class score_t, class sscore_t, size_t W>
void align_vec( persistent_state<score_t> &ps, size_t asize, const std::vector<uint8_t> &b, const scoring_matrix &m, aligned_buffer<sscore_t> &qprofile, const sscore_t gap_open, const sscore_t gap_extend, std::vector<int> &out ) {
 
    typedef vector_unit<score_t,W> vu;
 
    typedef typename vu::vec_t vec_t;
    
//     assert( a.size() % W == 0 );
//     const size_t asize = a.size() / W;
    
//     std::vector<score_t> s( a.size());
//     std::vector<score_t> si( a.size());
    
    aligned_buffer<score_t> &s = ps.s;
    aligned_buffer<score_t> &si = ps.si;
    
    aligned_buffer<score_t> pb(W);
    
    if( s.size() < asize * W ) {
        s.resize( asize * W );
        si.resize( asize * W );
    }
    const score_t bias = vu::BIAS;
    std::fill( s.begin(), s.end(), bias );
    std::fill( si.begin(), si.end(), bias );
    const score_t SMALL = vu::SMALL_VALUE;
    vec_t max_vec = vu::set1(SMALL);
    
    const vec_t GAP_EXT_vec = vu::set1(score_t(-gap_extend)); // these values are _subtracted_ from the score. so use negative of user parameters!
    const vec_t GAP_OPEN_vec = vu::set1(score_t(-gap_open));
            
    
//     vec_t len_vec = vu::load( len.m_ptr );
#define LOCAL_ALIGN
    for( size_t ib = 0; ib < b.size(); ib++ ) {
#ifndef LOCAL_ALIGN        
        const bool lastrow = ib == (b.size() - 1);
#endif        
        char bc = b[ib];
//         std::cout << "bc: " << bc << std::endl;
        assert( bc < char(qprofile.size() / (W * asize)));
        
        vec_t last_sl_vec = vu::set1(SMALL);
        vec_t last_sc_vec = vu::set1(vu::BIAS);
        vec_t last_sdiag_vec = vu::set1(vu::BIAS);
//         std::cout << "sbm: " << m.state_backmap(bc) << " " << (W * asize) << std::endl;
       // sscore_t *qpp_iter = qprofile(m.state_backmap(bc) * W * asize);
        sscore_t *qpp_iter = qprofile( bc * W * asize);
        score_t * __restrict s_iter = s.base();
        score_t * __restrict si_iter = si.base();
        score_t * __restrict s_end = s_iter + (asize * W);
        
        for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {
            
            const vec_t match_vec = vu::load( (score_t*) qpp_iter );
            
//             getchar();
            const vec_t sm_vec = vu::add( last_sdiag_vec, match_vec );
            
            last_sdiag_vec = vu::load( s_iter );
            
#ifdef LOCAL_ALIGN            
            const vec_t sm_zero_vec = vu::max( sm_vec, vu::set1(vu::BIAS));
#else
            const vec_t sm_zero_vec = sm_vec;
#endif
            
            const vec_t last_sc_OPEN_vec = vu::sub( last_sc_vec, GAP_OPEN_vec ); 
            
            const vec_t sl_vec = vu::max( vu::sub( last_sl_vec, GAP_EXT_vec ), last_sc_OPEN_vec );
//             vu::println( GAP_EXT_vec, pb.m_ptr ); getchar();
            
            
            last_sl_vec = sl_vec;
            

            const vec_t su_GAP_EXTEND_vec = vu::sub( vu::load(si_iter), GAP_EXT_vec );
            const vec_t su_vec = vu::max( su_GAP_EXTEND_vec, vu::sub( last_sdiag_vec, GAP_OPEN_vec) );
            
            vu::store( su_vec, si_iter );
            
            const vec_t sc_vec = vu::max( sm_zero_vec, vu::max( sl_vec, su_vec ));
            
//             vu::println( max_vec, pb.m_ptr );
            
            last_sc_vec = sc_vec;
            vu::store( sc_vec, s_iter ); 
#ifdef LOCAL_ALIGN  
            max_vec = vu::max( sc_vec, max_vec );
            
#else
            if( s_iter == s_end - W || lastrow ) {
                max_vec = vu::max( sc_vec, max_vec );
            }
#endif
        }
        
    }
    
    ps.out.resize(W);
    vu::store( max_vec, ps.out.base() );
    
    out.resize(0);
    out.reserve(W);
    for( size_t i = 0; i < W; i++ ) {
        out.push_back(int(ps.out[i]) - vu::BIAS );
//         std::cout << "x: " << int(ps.out.m_ptr[i]) << "\n";
    }
    //return max;
}


#endif