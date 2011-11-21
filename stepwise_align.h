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


#ifndef __stepwise_align_h
#define __stepwise_align_h
#include <iostream>
#include <ostream>
#include <iterator>
#include <cassert>

#include "aligned_buffer.h"
#include "vec_unit.h"
#include "parsimony.h"
#include "fasta.h"
//
// general: the freeshift version of this type of aligner differs from the
// implementation used in papara mainly by allowing free gaps also on the reference side
//
// update: the QS-side free-gaps seem to work quite well for normal papara as well, plus
// the positive/negaitve scoring scheme also works better than the original one from papara.
// Because of these factors this is bound to become the standard alignment kernel for 'papara 2.0'
//

template<typename score_t>
struct align_arrays {
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
};

template<typename score_t,bool global>
score_t align_pvec_score( std::vector<uint8_t> &a, std::vector<uint8_t> &a_aux, const std::vector<uint8_t> &b, score_t match_score, score_t match_cgap, score_t gap_open, score_t gap_extend, align_arrays<score_t> &arr ) {

    if( arr.s.size() < a.size()  ) {
        arr.s.resize( a.size() );
        arr.si.resize( a.size() );
    }
    
    
    if( global ) {
        for( size_t i = 0; i < a.size(); ++i ) {
            arr.s[i] = gap_open + (i+1) * gap_extend;
        }
    } else {
        std::fill( arr.s.begin(), arr.s.end(), 0 );
    }
    const score_t SMALL = -32000;
    std::fill( arr.si.begin(), arr.si.end(), SMALL );
    

    score_t max_score = SMALL;
    
    for( size_t ib = 0; ib < b.size(); ib++ ) {
        int bc = b[ib];
        
        score_t last_sl = SMALL;
        score_t last_sc;// = 0.0;
        score_t last_sdiag;// = 0.0;
        if( global ) {
            if( ib == 0 ) {
                last_sdiag = 0;
            } else {
                last_sdiag = gap_open + ib * gap_extend;
            }
            
            last_sc = gap_open + (ib+1) * gap_extend;
            
        } else {
            last_sc = 0;
            last_sdiag = 0;
        }
        
        
        score_t * __restrict s_iter = arr.s.base();
        score_t * __restrict si_iter = arr.si.base();
        score_t * __restrict s_end = s_iter + a.size();
        bool lastrow = ib == (b.size() - 1);
            
        for( size_t ia = 0; ia < a.size(); ++ia, ++s_iter, ++si_iter ) {  
            score_t ac = a[ia];
            const bool cgap = a_aux[ia] == AUX_CGAP;
            
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
            } else {
                sl = last_sc_OPEN;
            }
            
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
            
            if( global ) {
                max_score = sc;
            } else {
                if( s_iter == (s_end - 1) || lastrow ) {
                    if( sc  > max_score ) {
                        
                        max_score = sc;
                    }
                }
            }
        }
    }

    return max_score;
    
}



#if 0
template<typename score_t>
struct align_vec_arrays {
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
};

template<typename score_t, size_t W, bool global, typename bstate_t>
void align_pvec_score_vec( aligned_buffer<score_t> &a_prof, aligned_buffer<score_t> &a_aux_prof, const std::vector<bstate_t> &b, const score_t match_score_sc, const score_t match_cgap_sc, const score_t gap_open_sc, const score_t gap_extend_sc, aligned_buffer<score_t> &out, align_vec_arrays<score_t> &arr ) {


    typedef vector_unit<score_t,W> vu;
    typedef typename vu::vec_t vec_t;



    size_t av_size = a_prof.size();

//     if( arr.s.size() < av_size  ) {
    arr.s.resize( av_size );
    arr.si.resize( av_size );
//     }


    if( arr.s.size() % W != 0 ) {
        throw std::runtime_error( "profile length not multiple of vec-unit width." );
    }

    if( global ) {
        int i = 0;

        for( typename aligned_buffer<score_t>::iterator it = arr.s.begin(); it != arr.s.end(); it += W, ++i ) {
            std::fill( it, it + W, gap_open_sc + (i+1) * gap_extend_sc );
        }
    } else {
        std::fill( arr.s.begin(), arr.s.end(), 0 );
    }

    const score_t SMALL = vu::SMALL_VALUE;
    std::fill( arr.si.begin(), arr.si.end(), SMALL );


    vec_t max_score = vu::set1(SMALL);

    const vec_t zero = vu::setzero();
    const vec_t gap_extend = vu::set1(gap_extend_sc);
    const vec_t gap_open = vu::set1(gap_open_sc);
    const vec_t match_cgap = vu::set1( match_cgap_sc );
    const vec_t match_score = vu::set1( match_score_sc );




    for( size_t ib = 0; ib < b.size(); ++ib ) {
        vec_t bc = vu::set1(b[ib]);

        vec_t last_sl = vu::set1(SMALL);
        vec_t last_sc;// = vu::set1(0);
        vec_t last_sdiag; // = vu::set1(0);

        if( global ) {
            if( ib == 0 ) {
                last_sdiag = vu::set1(0);
            } else {
                last_sdiag = vu::set1(gap_open_sc + ib * gap_extend_sc);
            }

            last_sc = vu::set1(gap_open_sc + (ib+1) * gap_extend_sc);

        } else {
            last_sc = vu::set1(0);;
            last_sdiag = vu::set1(0);;
        }


        score_t * __restrict a_prof_iter = a_prof.base();
        score_t * __restrict a_aux_prof_iter = a_aux_prof.base();
        score_t * __restrict s_iter = arr.s.base();
        score_t * __restrict si_iter = arr.si.base();
        score_t * __restrict s_end = s_iter + av_size;
        bool lastrow = ib == (b.size() - 1);

        vec_t row_max_score = vu::set1(SMALL);

        for(; s_iter != s_end; a_prof_iter += W, a_aux_prof_iter += W, s_iter += W, si_iter += W ) {

            //break_aloop = s_iter == (s_end - W);


            const vec_t ac = vu::load( a_prof_iter );
            const vec_t cgap = vu::load( a_aux_prof_iter );

            const vec_t non_match = vu::cmp_eq( vu::bit_and( ac, bc ), zero );



            // match increase: sum of match and match_cgap score/penalty
            const vec_t sm_inc = vu::add( vu::bit_andnot(non_match, match_score), vu::bit_and(cgap, match_cgap) );

            // match score (last_sdiag is preloaded in the previous iteration)
            const vec_t sm = vu::add(last_sdiag, sm_inc );

            // last_sdiag now is preloaded with the 'above' value (it will become the 'diag'
            // element in the next iteration.)



            last_sdiag = vu::load( s_iter );

            // open score for gap-from-left (last_sc is the previous main-cell score)
            const vec_t sl_open = vu::add( last_sc, vu::bit_andnot( cgap, gap_open) );

            // extension score for gap-from-left (last_sl is the previous gap-from-left score)
            const vec_t sl_extend = vu::add( last_sl, vu::bit_andnot( cgap, gap_extend) );

            const vec_t sl = vu::max( sl_open, sl_extend );
            last_sl = sl;


            // open score for gap from above (remember, last_sdiag is preloaded with the 'above' value)
            const vec_t su_open = vu::add( last_sdiag, gap_open);

            // extension score for gap from above (si_iter is the previous gap-from-above score)
            const vec_t su_GAP_EXTEND = vu::add( vu::load( si_iter ), gap_extend );

            // chose between open/extension for gap-from-above
            const vec_t su = vu::max( su_GAP_EXTEND, su_open );// = max( su_GAP_EXTEND,  );
            vu::store( su, si_iter );


            // overall max = new main-cell score
            const vec_t sc = vu::max( sm, vu::max( su, sl ) );
            last_sc = sc;


            vu::store( sc, s_iter );

            if( !global ) {
                row_max_score = vu::max( row_max_score, sc );
            }
        }

        if( global ) {
            max_score = last_sc;
        } else {
            max_score = vu::max( max_score, last_sc );
            if( lastrow ) {
                max_score = vu::max( max_score, row_max_score );
            }
        }
    }


    vu::store( max_score, out.base() );
}
#else

template<typename score_t>
struct align_vec_arrays {
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;
};


template<typename score_t>
struct ali_ptr_block_t {
    score_t * __restrict a_prof_iter;
    score_t * __restrict a_aux_prof_iter;
    //score_t * __restrict s_iter;
    //score_t * __restrict si_iter;
};

template<typename vec_t>
struct ali_score_block_t {
    vec_t last_sl;
    vec_t last_sc;
    vec_t last_sdiag;
};

template<typename score_t, size_t W, bool global, typename bstate_t>
void align_pvec_score_vec( aligned_buffer<score_t> &a_prof, aligned_buffer<score_t> &a_aux_prof, const std::vector<bstate_t> &b, const score_t match_score_sc, const score_t match_cgap_sc, const score_t gap_open_sc, const score_t gap_extend_sc, aligned_buffer<score_t> &out, align_vec_arrays<score_t> &arr ) {
    
    assert( !global );
    
    typedef vector_unit<score_t,W> vu;
    typedef typename vu::vec_t vec_t;
    
    
        
    size_t av_size = a_prof.size();
    
//     if( arr.s.size() < av_size  ) {
    arr.s.resize( av_size );
    arr.si.resize( av_size );
//     }
    
    
    if( arr.s.size() % W != 0 ) {
        throw std::runtime_error( "profile length not multiple of vec-unit width." );
    }
    
    

    const score_t SMALL = vu::SMALL_VALUE;
    std::fill( arr.s.begin(), arr.s.end(), 0 );
    std::fill( arr.si.begin(), arr.si.end(), SMALL );
    

    vec_t max_score = vu::set1(SMALL);
    
    const vec_t zero = vu::setzero();
    const vec_t gap_extend = vu::set1(gap_extend_sc);
    const vec_t gap_open = vu::set1(gap_open_sc);
    const vec_t match_cgap = vu::set1( match_cgap_sc );
    const vec_t match_score = vu::set1( match_score_sc );

    score_t * const a_end_final = a_prof.base() + av_size;

    bool done = false;




    const size_t block_width = 512;

    assert( a_prof.size() >= block_width * W );

    ali_score_block_t<vec_t> btemp;

    btemp.last_sl = vu::set1(SMALL);
    btemp.last_sc = vu::set1(0);;
    btemp.last_sdiag = vu::set1(0);;

    std::vector<ali_score_block_t<vec_t> > blocks( b.size(), btemp );


    ali_ptr_block_t<score_t> ptr_block_outer;
    ptr_block_outer.a_prof_iter = a_prof.base();
    ptr_block_outer.a_aux_prof_iter = a_aux_prof.base();
    //ptr_block_outer.s_iter = arr.s.base();
    //ptr_block_outer.si_iter = arr.si.base();



    while( !done ) {
        

        ali_ptr_block_t<score_t> ptr_block = ptr_block_outer;

        std::fill( arr.s.begin(), arr.s.begin() + W * block_width, 0 );
        std::fill( arr.si.begin(), arr.si.begin() + W * block_width, SMALL );

        for( size_t ib = 0; ib < b.size(); ++ib ) {
            vec_t bc = vu::set1(b[ib]);



            
            
//            score_t * __restrict a_prof_iter = a_prof.base();
//            score_t * __restrict a_aux_prof_iter = a_aux_prof.base();
//            score_t * __restrict s_iter = arr.s.base();
//            score_t * __restrict si_iter = arr.si.base();

            bool lastrow = ib == (b.size() - 1);
            
            vec_t row_max_score = vu::set1(SMALL);

            ali_score_block_t<vec_t> block = blocks[ib];

            ptr_block = ptr_block_outer;
            score_t * a_end = ptr_block.a_prof_iter + W * block_width;

            if( a_end > a_end_final ) {
                a_end = a_end_final;
            }
            score_t * __restrict s_iter = arr.s.base();
            score_t * __restrict si_iter = arr.si.base();
            //ptr_block_outer.si_iter = arr.si.base();

            for(; ptr_block.a_prof_iter != a_end; ptr_block.a_prof_iter += W, ptr_block.a_aux_prof_iter += W, s_iter += W, si_iter += W ) {
                //break_aloop = s_iter == (s_end - W);
                const vec_t ac = vu::load( ptr_block.a_prof_iter );
                const vec_t cgap = vu::load( ptr_block.a_aux_prof_iter );

                const vec_t non_match = vu::cmp_eq( vu::bit_and( ac, bc ), zero );



                // match increase: sum of match and match_cgap score/penalty
                const vec_t sm_inc = vu::add( vu::bit_andnot(non_match, match_score), vu::bit_and(cgap, match_cgap) );

                // match score (last_sdiag is preloaded in the previous iteration)
                const vec_t sm = vu::add(block.last_sdiag, sm_inc );

                // last_sdiag now is preloaded with the 'above' value (it will become the 'diag'
                // element in the next iteration.)



                block.last_sdiag = vu::load( s_iter );

                // open score for gap-from-left (last_sc is the previous main-cell score)
                const vec_t sl_open = vu::add( block.last_sc, vu::bit_andnot( cgap, gap_open) );

                // extension score for gap-from-left (last_sl is the previous gap-from-left score)
                const vec_t sl_extend = vu::add( block.last_sl, vu::bit_andnot( cgap, gap_extend) );

                const vec_t sl = vu::max( sl_open, sl_extend );
                block.last_sl = sl;


                // open score for gap from above (remember, last_sdiag is preloaded with the 'above' value)
                const vec_t su_open = vu::add( block.last_sdiag, gap_open);

                // extension score for gap from above (si_iter is the previous gap-from-above score)
                const vec_t su_GAP_EXTEND = vu::add( vu::load( si_iter ), gap_extend );

                // chose between open/extension for gap-from-above
                const vec_t su = vu::max( su_GAP_EXTEND, su_open );// = max( su_GAP_EXTEND,  );
                vu::store( su, si_iter );


                // overall max = new main-cell score
                const vec_t sc = vu::max( sm, vu::max( su, sl ) );
                block.last_sc = sc;


                vu::store( sc, s_iter );

                if( !global ) {
                    row_max_score = vu::max( row_max_score, sc );
                }
            }

            done = ptr_block.a_prof_iter == a_end_final;

            if( done ) {
                max_score = vu::max( max_score, block.last_sc );
            }
            if( lastrow ) {
                max_score = vu::max( max_score, row_max_score );
            }

            blocks[ib] = block;

        }

        ptr_block_outer = ptr_block;
    }
    
    vu::store( max_score, out.base() );
}
#endif

//
// the 'full enchilada' aligner including traceback.
// The freeshift/global implementations are currently separate mainly because they
// act as reference for the optimized versions and I don't want to screw them up too much...
//

template<typename score_t>
struct align_arrays_traceback {
    aligned_buffer<score_t> s;
    aligned_buffer<score_t> si;  
    
    std::vector<uint8_t> tb;
};


template<typename score_t, typename aiter, typename auxiter, typename biter>
score_t align_freeshift_pvec( aiter astart, aiter aend, auxiter auxstart, biter bstart, biter bend, score_t match_score, score_t match_cgap, score_t gap_open, score_t gap_extend, std::vector<uint8_t>& tb_out, align_arrays_traceback<score_t> &arr ) {

    const uint8_t b_sl_stay = 0x1;
    const uint8_t b_su_stay = 0x2;
    const uint8_t b_s_l = 0x4;
    const uint8_t b_s_u = 0x8;


    const size_t asize = std::distance(astart, aend);
    const size_t bsize = std::distance(bstart, bend);

    if( arr.s.size() < asize  ) {
        arr.s.resize( asize );
        arr.si.resize( asize );
    }

    std::fill( arr.s.begin(), arr.s.end(), 0 );
    std::fill( arr.si.begin(), arr.si.end(), 0 );
    const score_t SMALL = -32000;

    score_t max_score = SMALL;
    int max_a = 0;
    int max_b = 0;

    arr.tb.resize( asize * bsize );


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

    index_calc ic( asize, bsize );


    for( size_t ib = 0; ib < bsize; ib++ ) {
        int bc = *(bstart + ib);

        score_t last_sl = SMALL;
        score_t last_sc = 0.0;
        score_t last_sdiag = 0.0;

        score_t * __restrict s_iter = arr.s.base();
        score_t * __restrict si_iter = arr.si.base();
        score_t * __restrict s_end = s_iter + asize;
        bool lastrow = ib == (bsize - 1);

        //for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {

        for( size_t ia = 0; ia < asize; ++ia, ++s_iter, ++si_iter ) {
            //score_t match = sm.get_score( a[ia], bc );
            uint8_t tb_val = 0;
            int ac = *(astart + ia);
            const bool cgap = *(auxstart + ia)  == AUX_CGAP;

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
                //arr.sl_stay[ic(ia,ib)] = true;
                tb_val |= b_sl_stay;
            } else {
                sl = last_sc_OPEN;
            }

            last_sl = sl;


            score_t su_gap_open = last_sdiag + gap_open;
            score_t su_GAP_EXTEND = *si_iter + gap_extend;

            score_t su;// = max( su_GAP_EXTEND,  );
            if( su_GAP_EXTEND > su_gap_open ) {
                su = su_GAP_EXTEND;
                //arr.su_stay[ic(ia,ib)] = true;
                tb_val |= b_su_stay;
            } else {
                su = su_gap_open;
            }


            *si_iter = su;

            score_t sc;
            if( (su > sl) && su > sm ) {
                sc = su;
                //arr.s_u[ic(ia,ib)] = true;
                tb_val |= b_s_u;
            } else if( ( sl >= su ) && sl > sm ) {
                sc = sl;
                //arr.s_l[ic(ia,ib)] = true;
                tb_val |= b_s_l;
            } else { // implicit: sm_zero > sl && sm_zero > su
                sc = sm;
            }


            last_sc = sc;
            *s_iter = sc;
            arr.tb[ic(ia,ib)] = tb_val;

            if( s_iter == s_end - 1 || lastrow ) {
                if( sc > max_score ) {
                    max_a = int(ia);
                    max_b = int(ib);
                    max_score = sc;
                }


            }
        }

    }

//    std::cout << "max score: " << max_score << "\n";
//    std::cout << "max " << max_a << " " << max_b << "\n";

    int ia = asize - 1;
    int ib = bsize - 1;

    assert( ia == max_a || ib == max_b );

    bool in_l = false;
    bool in_u = false;

    while( ia > max_a ) {
        tb_out.push_back(1);
        --ia;
    }

    while( ib > max_b ) {
        tb_out.push_back(2);
        --ib;
    }

    while( ia >= 0 && ib >= 0 ) {
        size_t c = ic( ia, ib );

        if( !in_l && !in_u ) {
            in_l = (arr.tb[c] & b_s_l) != 0;
            in_u = (arr.tb[c] & b_s_u) != 0;

            if( !in_l && !in_u ) {
                tb_out.push_back(0);
                --ia;
                --ib;
            }

        }

        if( in_u ) {
            tb_out.push_back(2);
            --ib;

            in_u = (arr.tb[c] & b_su_stay) != 0;
        } else if( in_l ) {
            tb_out.push_back(1);
            --ia;

            in_l = (arr.tb[c] & b_sl_stay) != 0;
        }


    }

    while( ia >= 0 ) {
        tb_out.push_back(1);
        --ia;
    }

    while( ib >= 0 ) {
        tb_out.push_back(2);
        --ib;
    }
    return max_score;
}

// backward compatibility wrapper
template<typename score_t, typename state_t>
inline score_t align_freeshift_pvec( std::vector<state_t> &a, std::vector<state_t> &a_aux, std::vector<state_t> &b, score_t match_score, score_t match_cgap, score_t gap_open, score_t gap_extend, std::vector<uint8_t>& tb_out, align_arrays_traceback<score_t> &arr ) {
    return align_freeshift_pvec( a.begin(), a.end(), a_aux.begin(), b.begin(), b.end(), match_score, match_cgap, gap_open, gap_extend, tb_out, arr );
}



template<typename score_t>
score_t align_global_pvec( std::vector<uint8_t> &a, std::vector<uint8_t> &a_aux, std::vector<uint8_t> &b, score_t match_score, score_t match_cgap, score_t gap_open, score_t gap_extend, std::vector<uint8_t>& tb_out, align_arrays_traceback<score_t> &arr ) {
  
    const uint8_t b_sl_stay = 0x1;
    const uint8_t b_su_stay = 0x2;
    const uint8_t b_s_l = 0x4;
    const uint8_t b_s_u = 0x8;
    
    
    if( arr.s.size() < a.size()  ) {
        arr.s.resize( a.size() );
        arr.si.resize( a.size() );
    }
    
    
    const score_t SMALL = -32000;
    //fill( arr.s.begin(), arr.s.end(), 0 );
    for( size_t i = 0; i < a.size(); ++i ) {
        arr.s[i] = gap_open + (i+1) * gap_extend;
    }
    
    std::fill( arr.si.begin(), arr.si.end(), SMALL );
    

    score_t max_score = SMALL;
    
    
    arr.tb.resize( a.size() * b.size() );
    
    
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
        
        score_t last_sdiag;
        score_t last_sc = gap_open + (ib+1) * gap_extend;
        if( ib == 0 ) {
            last_sdiag = 0;
        } else {
            last_sdiag = gap_open + ib * gap_extend;
        }
        
        score_t * __restrict s_iter = arr.s.base();
        score_t * __restrict si_iter = arr.si.base();
        
        
        
        //for( ; s_iter != s_end; s_iter += W, si_iter += W, qpp_iter += W ) {
            
        for( size_t ia = 0; ia < a.size(); ++ia, ++s_iter, ++si_iter ) {  
            //score_t match = sm.get_score( a[ia], bc );
            uint8_t tb_val = 0;
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
                //arr.sl_stay[ic(ia,ib)] = true;
                tb_val |= b_sl_stay;
            } else {
                sl = last_sc_OPEN;
            }
            
            last_sl = sl;
            

            score_t su_gap_open = last_sdiag + gap_open;
            score_t su_GAP_EXTEND = *si_iter + gap_extend;
            
            score_t su;// = max( su_GAP_EXTEND,  );
            if( su_GAP_EXTEND > su_gap_open ) {
                su = su_GAP_EXTEND;
                //arr.su_stay[ic(ia,ib)] = true;
                tb_val |= b_su_stay;
            } else {
                su = su_gap_open;
            }
            
            
            *si_iter = su;
            
            score_t sc;
            if( (su > sl) && su > sm ) {
                sc = su;
                //arr.s_u[ic(ia,ib)] = true;
                tb_val |= b_s_u;
            } else if( ( sl >= su ) && sl > sm ) {
                sc = sl;
                //arr.s_l[ic(ia,ib)] = true;
                tb_val |= b_s_l;
            } else { // implicit: sm_zero > sl && sm_zero > su
                sc = sm;
            }

            
            last_sc = sc;
            *s_iter = sc;
            arr.tb[ic(ia,ib)] = tb_val;
            
            
            max_score = sc;

        }
        
    }
    
//     std::cout << "max score: " << max_score << "\n";
    
    
    int ia = a.size() - 1;
    int ib = b.size() - 1;

    
    bool in_l = false;
    bool in_u = false;
    
//     while( ia > max_a ) {
//         tb_out.push_back(1);
//         --ia;
//     }
//     
//     while( ib > max_b ) {
//         tb_out.push_back(2);
//         --ib;
//     }
    
    while( ia >= 0 && ib >= 0 ) {
        size_t c = ic( ia, ib );
        
        if( !in_l && !in_u ) {
            in_l = (arr.tb[c] & b_s_l) != 0;
            in_u = (arr.tb[c] & b_s_u) != 0;
            
            if( !in_l && !in_u ) {
                tb_out.push_back(0);
                --ia;
                --ib;
            }
            
        }
        
        if( in_u ) {
            tb_out.push_back(2);
            --ib;
            
            in_u = (arr.tb[c] & b_su_stay) != 0;
        } else if( in_l ) {
            tb_out.push_back(1);
            --ia;
            
            in_l = (arr.tb[c] & b_sl_stay) != 0;
        }
        
        
    }
    
    while( ia >= 0 ) {
        tb_out.push_back(1);
        --ia;
    }
    
    while( ib >= 0 ) {
        tb_out.push_back(2);
        --ib;
    }
    return max_score;
}


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
    int max_a = 0;
    int max_b = 0;
    
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
                    max_a = int(ia);
                    max_b = int(ib);
                    max_score = sc;
                }
                
                
            }
        }
        
    }
    
//    std::cout << "max score: " << max_score << "\n";
//    std::cout << "max " << max_a << " " << max_b << "\n";
//    std::cout << "size " << a.size() << " " << b.size() << "\n";
    
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
    }
    
    while( ib > max_b ) {
        a_tb.push_back('-');
        b_tb.push_back(b[ib]);
        --ib;
    }
    
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
    
//     std::copy( a.begin(), a.end(), std::ostream_iterator<char>(std::cout) );
//     std::cout << "\n";
//     std::copy( b.begin(), b.end(), std::ostream_iterator<char>(std::cout) );
//     std::cout << "\n";
}

#endif
