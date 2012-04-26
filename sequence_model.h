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

#ifndef SEQUENCE_MODEL_H_
#define SEQUENCE_MODEL_H_

#include <stdint.h>
#include <cctype>
#include <cstddef>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "ivymike/algorithm.h"

#ifndef _MSC_VER
template<typename T>
size_t popcount( T v ) {
	return __builtin_popcount(v);
}
#else
#include <intrin.h>

inline size_t popcount( unsigned short v ) {
	return __popcnt16(v);
}

inline size_t popcount( unsigned int v ) {
	return __popcnt(v);
}

inline size_t popcount( unsigned __int64 v ) {
	return __popcnt64(v);
}


#endif


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

    const static std::vector<char> inverse_meaning;
//    const static std::vector<pars_state_t> bit_vector;

    static uint8_t normalize( uint8_t c ) {
        c = std::toupper(c);

        switch( c ) {
        case 'U':
            return 'T';
        case 'N':
            return '-';
        default:
            return c;
        }

    }


    static pars_state_t s2p( uint8_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            std::cerr << "illegal character: " << int(c) << "\n";
            throw std::runtime_error( "illegal character in DNA/RNA sequence");
        }

        return pars_state_t(idx);
    }

    static uint8_t s2c( uint8_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            std::cerr << "illegal character: " << int(c) << "\n";
            throw std::runtime_error( "illegal character in DNA/RNA sequence");
        }

        return idx;
    }

    static pars_state_t c2p( uint8_t c ) {
        return pars_state_t(c);
    }


    static uint8_t p2s( pars_state_t c ) {
        return inverse_meaning.at(c);
    }


    static inline bool is_single(pars_state_t ps) {
        return !is_gap(ps) && ps != 0;
    }

    static inline bool is_gap(pars_state_t ps) {
        return ps == gap_state();
    }

    static inline bool cstate_is_gap( uint8_t cs) {
       return cs == inverse_meaning.size() - 1;
   }

    static inline pars_state_t gap_state() {
        return inverse_meaning.size() - 1;
    }

    static inline size_t num_cstates() {
        return inverse_meaning.size();
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

    static uint8_t s2c( uint8_t c ) {
        c = normalize(c);
        ptrdiff_t idx = std::distance(inverse_meaning.begin(),
                                   std::find(inverse_meaning.begin(), inverse_meaning.end(), c ) );

        assert( idx >= 0 );

        if( size_t(idx) >= inverse_meaning.size() ) {
            std::cerr << "illegal character: " << int(c) << "\n";
            throw std::runtime_error( "illegal character in DNA/RNA sequence");
        }

        return idx;
    }

    static pars_state_t c2p( uint8_t c ) {
        return bit_vector.at(c);
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
        return popcount(ps) == 1;
    }

    static inline bool is_gap(pars_state_t ps) {
        return ps == gap_state();
    }

    static inline bool cstate_is_gap( uint8_t cs) {
       return cs == inverse_meaning.size() - 1;
   }


    static inline pars_state_t gap_state() {
        return bit_vector.back(); // by convention the last element of bit_vector is the gap state
    }

    static inline size_t num_cstates() {
        return inverse_meaning.size();
    }
};


}



#endif /* SEQUENCE_MODEL_H_ */
