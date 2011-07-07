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

#ifndef __pairwise_seq_distance_h
#define __pairwise_seq_distance_h


// under certain circumstances gcc produces more consistent performance when the heavy-lifting is done in included code. maybe due to better optimization inside a compilation unit.
#ifdef PWDIST_INLINE
#define PSD_DECLARE_INLINE inline
#include "pairwise_seq_distance.cpp"
#else
#include <cstddef>
#include <vector>
#include <stdint.h>
#include "ivymike/tdmatrix.h"

class scoring_matrix;
bool pairwise_seq_distance( const std::vector< std::vector<uint8_t> > &seq_raw1, const std::vector< std::vector<uint8_t> > &seq_raw2, bool identical, ivy_mike::tdmatrix<int> &out_scores, scoring_matrix &sm, const int gap_open, const int gap_extend, const size_t n_thread );

#endif
inline bool pairwise_seq_distance( const std::vector< std::vector<uint8_t> > &seq_raw, ivy_mike::tdmatrix<int> &out_scores, scoring_matrix &sm, const int gap_open, const int gap_extend, const size_t n_thread ) {
    return pairwise_seq_distance( seq_raw, seq_raw, true, out_scores, sm, gap_open, gap_extend, n_thread );
}

#endif