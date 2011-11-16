/*
 * sequence_model.cpp
 *
 *  Created on: 11/11/2011
 *      Author: sim
 */


#include "sequence_model.h"


using sequence_model::model;
using sequence_model::tag_aa;
using sequence_model::tag_dna;

namespace raxml_aa_meaning {
// is it officially legal to initialize const static members in the header? I guess c++ removes redundant definitions during linking...
const static char inverse[23] = {'A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', '-'};
const static unsigned int bitVector[23] = {1, 2, 4, 8, 16, 32, 64, 128,
                                      256, 512, 1024, 2048, 4096,
                                      8192, 16384, 32768, 65536, 131072, 262144,
                                      524288, 12 /* N | D */, 96 /*Q | E*/, 1048575 /* - */};
}

const std::vector<char> model<tag_aa>::inverse_meaning(raxml_aa_meaning::inverse, raxml_aa_meaning::inverse + ivy_mike::arrlen(raxml_aa_meaning::inverse));
const std::vector<unsigned int> model<tag_aa>::bit_vector(raxml_aa_meaning::bitVector, raxml_aa_meaning::bitVector + ivy_mike::arrlen(raxml_aa_meaning::bitVector));



namespace raxml_dna_meaning {
const char inverse[16]   = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'};
}



