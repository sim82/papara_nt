/*
 * sequence_model.h
 *
 *  Created on: 10/11/2011
 *      Author: sim
 */

#ifndef SEQUENCE_MODEL_H_
#define SEQUENCE_MODEL_H_

#include <stdint.h>
#include <cctype>
#include <cstddef>
#include <cassert>
#include <iostream>
#include "ivymike/algorithm.h"

namespace sequence_model {


class tag_dna;
class tag_aa;


template<typename TAG>
class model {
    //static uint8_t normalize( uint8_t c );
};

template<>
class model<tag_dna> {
public:

    typedef uint8_t pars_state_t;


    static uint8_t normalize( uint8_t c ) {
        c = std::toupper(c);

        if( c == 'U' ) {
            c = 'T';
        }


        return c;
    }


    static pars_state_t s2p( uint8_t c ) {

        switch( c ) {
        case 'A':
        case 'a':
            return 0x1;

        case 'C':
        case 'c':
            return 0x2;

        case 'G':
        case 'g':
            return 0x4;

        case 'U':
        case 'u':
        case 'T':
        case 't':
            return 0x8;

        default:
            break;
        };
        return 0xf;
    }

    static uint8_t p2s( pars_state_t c ) {
        switch( c ) {
        case 0x1:
            return 'A';
        case 0x2:
            return 'C';
        case 0x4:
            return 'G';
        case 0x8:
            return 'T';
        case 0xf:
            return '-';

        default:
            break;

        };
        return 'X';
    }


    static inline bool is_single(pars_state_t ps) {
        return ps == 0x1 || ps == 0x2 || ps == 0x4 || ps == 0x8;
    }

    static inline bool is_gap(pars_state_t ps) {
        return ps == 0xf;
    }

    static inline pars_state_t gap_state() {
        return 0xf;
    }

};


// FIXME: fake aa model
template<>
class model<tag_aa> {
public:


    typedef uint32_t pars_state_t;
//    const static char inverseMeaningPROT[23];
//    const static unsigned int bitVectorAA[23];

    const static std::vector<char> inverse_meaning;
    const static std::vector<unsigned int> bit_vector;

    static inline uint8_t normalize( uint8_t c ) {
        return std::toupper(c);
    }


    static pars_state_t s2p( uint8_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );


//        std::cout << idx << "\n";

        // TODO: is there any reason to use more verbose range checking than '.at'?
        return bit_vector.at(idx);

    }

    static uint8_t p2s( pars_state_t c ) {

        ptrdiff_t idx = std::distance(bit_vector.begin(),
                                           std::find(bit_vector.begin(), bit_vector.end(), c ) );



        assert( idx >= 0 );
        if( idx >= ptrdiff_t(inverse_meaning.size()) ) {
            return 'X'; // parsimony state not representable as sequence character.
        }

        return inverse_meaning.at(idx);
    }


    static inline bool is_single(pars_state_t ps) {
        return __builtin_popcount(ps) == 1;
    }

    static inline bool is_gap(pars_state_t ps) {
        return ps == gap_state();
    }

    static inline pars_state_t gap_state() {
        return bit_vector.back(); // by convention the last element of bit_vector is the gap state
    }
};


}



#endif /* SEQUENCE_MODEL_H_ */
