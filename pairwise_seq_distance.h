#ifndef __pairwise_seq_distance_h
#define __pairwise_seq_distance_h



#ifdef PWDIST_INLINE
#define PSD_DECLARE_INLINE inline
#include "pairwise_seq_distance.cpp"
#else
#include <cstddef>
#include <vector>
#include <stdint.h>
#include "ivymike/tdmatrix.h"

class scoring_matrix;
void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq_raw, ivy_mike::tdmatrix<int> &out_scores, scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );
#endif


#endif