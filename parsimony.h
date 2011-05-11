#ifndef __parsimony_h
#define __parsimony_h

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <cassert>

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
void visit_edges( typename visitor::lnode *n, visitor &v, bool go_back = true ) {
    assert( n->back != 0 );
    
    v( n, n->back );
    
    if( go_back && n->back != 0 ) {
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


#endif