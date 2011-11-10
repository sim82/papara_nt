/*
 * sequence_model.h
 *
 *  Created on: 10/11/2011
 *      Author: sim
 */

#ifndef SEQUENCE_MODEL_H_
#define SEQUENCE_MODEL_H_



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
    typedef int16_t vu_scalar_t;
    const static size_t vu_width = 8;

    typedef uint8_t pars_state_t;


    static uint8_t normalize( uint8_t c ) {
        c = toupper(c);

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

};


// FIXME: fake aa model
template<>
class model<tag_aa> {
public:
    typedef int32_t vu_scalar_t;
    const static size_t vu_width = 4;

    typedef uint8_t pars_state_t;


    static uint8_t normalize( uint8_t c ) {
        c = toupper(c);

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

};



}



#endif /* SEQUENCE_MODEL_H_ */
