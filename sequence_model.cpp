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


#include "sequence_model.h"


using sequence_model::model;
using sequence_model::tag_aa;
using sequence_model::tag_dna;

namespace raxml_aa_meaning {
// is it officially legal to initialize const static members in the header? I guess c++ removes redundant definitions during linking...
const static char inverse[24] = {'A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', '-', 'X'};
const static unsigned int bitVector[24] = {1, 2, 4, 8, 16, 32, 64, 128,
                                      256, 512, 1024, 2048, 4096,
                                      8192, 16384, 32768, 65536, 131072, 262144,
                                      524288, 12 /* N | D */, 96 /*Q | E*/, 1048575 /* - */, 1048575 /* X */};
}

const std::vector<char> model<tag_aa>::inverse_meaning(raxml_aa_meaning::inverse, raxml_aa_meaning::inverse + ivy_mike::arrlen(raxml_aa_meaning::inverse));
const std::vector<unsigned int> model<tag_aa>::bit_vector(raxml_aa_meaning::bitVector, raxml_aa_meaning::bitVector + ivy_mike::arrlen(raxml_aa_meaning::bitVector));



namespace raxml_dna_meaning {
const char inverse[16]   = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-' };
//const char bitvector[17]   = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15 };
}



const std::vector<char> model<tag_dna>::inverse_meaning(raxml_dna_meaning::inverse, raxml_dna_meaning::inverse + ivy_mike::arrlen(raxml_dna_meaning::inverse));
//const std::vector<uint8_t> model<tag_dna>::bit_vector(raxml_dna_meaning::bitvector, raxml_dna_meaning::bitvector + ivy_mike::arrlen(raxml_dna_meaning::bitvector));
