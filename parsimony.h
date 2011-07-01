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

#ifndef __parsimony_h
#define __parsimony_h

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <cassert>
#include <stdint.h>
#if WIN32
#include <boost/shared_ptr.hpp>

namespace sptr = boost;
#else
#include <tr1/memory>
namespace sptr = std::tr1;
#endif

enum tip_case {
    TIP_TIP,
    TIP_INNER,
    INNER_INNER
    
};

enum aux_data {
    AUX_CGAP = 0x1,
    AUX_OPEN = 0x2
};

template<class lnode>
struct rooted_bifurcation {
    
    
    lnode * parent;
    lnode * child1;
    lnode * child2;
    
        
    tip_case tc;
    
    rooted_bifurcation() : parent(0), child1(0), child2(0) {}
    
    rooted_bifurcation( lnode *p, lnode *c1, lnode *c2, tip_case t ) 
        : parent(p),
        child1(c1),
        child2(c2),
        tc(t)
    {}
};


template<class lnode>
inline std::ostream &operator<<( std::ostream &os, const rooted_bifurcation<lnode> &rb ) {
    const char *tc;
    
    switch( rb.tc ) {
    case TIP_TIP:
        tc = "TIP_TIP";
        break;
        
    case TIP_INNER:
        tc = "TIP_INNER";
        break;
        
        case INNER_INNER:
        tc = "INNER_INNER";
        break;
    }
    
    os << tc << " " << *(rb.parent->m_data) << " " << *(rb.child1->m_data) << " " << *(rb.child2->m_data);
    return os;
}

template <class lnode, class container>
void rooted_traveral_order_rec( lnode *n, container &cont, bool incremental = false ) {
    
    
    
    lnode *n1 = n->next->back;
    lnode *n2 = n->next->next->back;
    
    assert( n->m_data->isTip || n1 != 0 );
    assert( n->m_data->isTip || n2 != 0 );
    
    n->towards_root = true;
    n->next->towards_root = false;
    n->next->next->towards_root = false;
    
    
    
    if( n1->m_data->isTip && n2->m_data->isTip ) {
        cont.push_front( rooted_bifurcation<lnode>( n, n1, n2, TIP_TIP ));
    } else if( n1->m_data->isTip && !n2->m_data->isTip ) {
        cont.push_front( rooted_bifurcation<lnode>( n, n1, n2, TIP_INNER ));
        
        if( !incremental || !n2->towards_root ) {
            rooted_traveral_order_rec( n2, cont );
        }
    } else if( !n1->m_data->isTip && n2->m_data->isTip ) {
        cont.push_front( rooted_bifurcation<lnode>( n, n2, n1, TIP_INNER ));
        
        if( !incremental || !n1->towards_root ) {
            rooted_traveral_order_rec( n1, cont );    
        }
        
    } else {
        cont.push_front( rooted_bifurcation<lnode>( n, n1, n2, INNER_INNER ));
        
        if( !incremental || !n1->towards_root ) {
            rooted_traveral_order_rec( n1, cont );
        }
        if( !incremental || !n2->towards_root ) {
            rooted_traveral_order_rec( n2, cont );
        }
    }
}

template <class lnode, class container>
void rooted_traveral_order( lnode *n1, lnode *n2, container &cont, bool incremental ) {
    
    if( !n1->m_data->isTip ) {
        rooted_traveral_order_rec<lnode, container>( n1, cont, incremental );
    }
    if( !n2->m_data->isTip ) {
        rooted_traveral_order_rec<lnode, container>( n2, cont, incremental );
    }
    
    
    //std::reverse( cont.begin(), cont.end());
}

template <class lnode>
lnode *towards_tree( lnode *n ) {
    int ct = 0;
    
    while( n->back == 0 ) {
        n = n->next;
        
        if( ct > 3 ) {
        
            throw std::runtime_error( "node not connected to tree" );
        }
        
        ct++;
    }
    
    return n;
    
}

template <class visitor>
void visit_lnode( typename visitor::lnode *n, visitor &v, bool go_back = true ) {
    v( n );
    
    if( go_back && n->back != 0 ) {
        visit_lnode( n->back, v, false );
    }
    if( n->next->back != 0 ) {
        visit_lnode( n->next->back, v, false );   
    }

    if( n->next->next->back != 0 ) {
        visit_lnode( n->next->next->back, v, false );
    }
};

template <class LNODE, class CONT = std::vector<sptr::shared_ptr<LNODE> > >
struct tip_collector {
	typedef LNODE lnode;
	typedef CONT container;
	
    //container<lnode *> m_nodes;
  
	container m_nodes;
	
public:
    void operator()( lnode *n ) {
        if( n->m_data->isTip ) {
            m_nodes.push_back(n->get_smart_ptr().lock());
        }
    }
};


template <class visitor>
void visit_edges( typename visitor::lnode *n, visitor &v, bool at_root = true ) {
    assert( n->back != 0 );
    
    
    // at the root, the edge between n and n->back will be visited when recursing to n->back
    if( !at_root ) {
        v( n, n->back );
    }
    
    if( at_root && n->back != 0 ) {
        visit_edges( n->back, v, false );
    }
    if( n->next->back != 0 ) {
        visit_edges( n->next->back, v, false );   
    }

    if( n->next->next->back != 0 ) {
        visit_edges( n->next->next->back, v, false );
    }
};

template <class LNODE>
struct edge_collector {
    typedef LNODE lnode;
    
    typedef std::pair<LNODE *, LNODE *> edge;
    std::vector<edge> m_edges;
    
public:
    void operator()( lnode *n1, lnode *n2 ) {
//         std::cout << "edge: " << n1 << " " << n2 << "\n";
        m_edges.push_back( edge( n1, n2 ) );
        
    }
};








//
// parsimony related stuff
//

typedef int parsimony_state;


// this class is so strange, we could name a new design-pattern after it.
// it does something like a static {} block in java...
class dna_parsimony_mapping_real {
    std::vector<uint8_t> m_p2d;
    std::vector<parsimony_state> m_d2p;
	char my_tolower( char x ) {
		// std::tolower needs std::locale. I will not touch i18n

		if( x >= 'A' && x <= 'Z' ) {
			return x - ('A' - 'a');
		} else {
			return x;
		}
		
	}
public:
    dna_parsimony_mapping_real() : m_p2d(16), m_d2p(256, -1)
    {
        const uint8_t  pd[18] =  {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', '-', 'N'};
        const uint32_t pds[18] = { 0,   1,   2,   3,   4,   5,   6,   7,   8,    8,   9,  10,  11,  12,  13,  14,  15,  15};
        //m_p2d.assign( pd, pd + 16 );


        // this is some weird code, which is supposed tp setup the dna->pstate and pstate->dna maps from ps and pds...
        for( int i = 0; i < 18; i++ ) {
            if( pds[i] >= m_p2d.size() ) {
                throw std::runtime_error( "dna/parsimony map initialization fsck'ed up" );
            }

            m_d2p[pd[i]] = pds[i];
            m_d2p[my_tolower(pd[i])] = pds[i];

            if( m_p2d[pds[i]] == 0 ) {
                m_p2d[pds[i]] = pd[i];
            }
        }

    }

    static dna_parsimony_mapping_real s_pdm;

    static uint8_t p2d( parsimony_state c ) {
        if( c < 0 || c > 15 ) {
            throw std::runtime_error( "illegal parsimony state" );
        }

        return s_pdm.m_p2d[c];
    }

    static parsimony_state d2p( uint8_t c ) {
        if( s_pdm.m_d2p[c] == -1 ) {

//             std::cerr << "illegal: " << c << "\n";
//             throw std::runtime_error( "illegal dna character" );
            return 0xf; // default to undefined pars state
        }

        return s_pdm.m_d2p[c];
    }

    static int d2aux( uint8_t c ) {
        parsimony_state ps = d2p(c);
        if( ps == 0xf ) {
            return AUX_CGAP;
        } else {
            return 0;
        }
    }

    static bool is_gap( uint8_t c ) {
        return d2p(c) == 0xf;
    }
};

dna_parsimony_mapping_real dna_parsimony_mapping_real::s_pdm;


struct dna_parsimony_mapping_simple {
    static parsimony_state d2p( uint8_t c ) {

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
            return 0xf;
        };
    }

    static uint8_t p2d( parsimony_state c ) {
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
            return 'X';
        };
    }
    static int d2aux( uint8_t c ) {
        parsimony_state ps = d2p(c);
        if( ps == 0xf ) {
            return AUX_CGAP;
        } else {
            return 0;
        }
    }
};

typedef dna_parsimony_mapping_real dna_parsimony_mapping;


#endif