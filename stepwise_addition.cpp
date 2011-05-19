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

#include "pairwise_seq_distance.h"

#include "ivymike/tdmatrix.h"


#include "fasta.h"



void pairwise_seq_distance(std::vector< std::vector<uint8_t> > &seq_raw, ivy_mike::tdmatrix<int> &, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );

int main() {
    //mapped_file qsf( "test_1604/1604.fa" );
	std::ifstream qsf( "test_1604/1604.fa" );
    std::vector<std::string> qs_names;
    std::vector<std::vector<uint8_t> > qs_seqs;
    
    std::vector<std::vector <uint8_t> > qs_nongappy;
    
    
    
    
    read_fasta( qsf, qs_names, qs_seqs);
    scoring_matrix sm( 3, 0 );
    
    ivy_mike::tdmatrix<int> out_scores(qs_seqs.size(), qs_seqs.size());
    pairwise_seq_distance(qs_seqs, out_scores, sm, -5, -2, 4);
    
}