#ifndef __parsimony_h
#define __parsimony_h

#include <algorithm>
#include <stdexcept>

template<class lnode>
struct rooted_bifurcation {
    enum tip_case {
        TIP_TIP,
        TIP_INNER,
        INNER_INNER
        
    };
    
    lnode * parent;
    lnode * child1;
    lnode * child2;
    
    tip_case tc;
    
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
    case rooted_bifurcation<lnode>::TIP_TIP:
        tc = "TIP_TIP";
        break;
        
    case rooted_bifurcation<lnode>::TIP_INNER:
        tc = "TIP_INNER";
        break;
        
        case rooted_bifurcation<lnode>::INNER_INNER:
        tc = "INNER_INNER";
        break;
    }
    
    os << tc << " " << *rb.parent << " " << *rb.child1 << " " << *rb.child2;
    
}

template <class lnode, class container>
void rooted_traveral_order_rec( lnode *n, container &cont ) {
    lnode *n1 = n->next->back;
    lnode *n2 = n->next->next->back;
    
    if( n1->m_data->isTip && n2->m_data->isTip ) {
        cont.push_back( rooted_bifurcation<lnode>( n, n1, n2, rooted_bifurcation<lnode>::TIP_TIP ));

    } else if( n1->m_data->isTip && !n2->m_data->isTip ) {
        cont.push_back( rooted_bifurcation<lnode>( n, n1, n2, rooted_bifurcation<lnode>::TIP_INNER ));
        rooted_traveral_order_rec( n2, cont );
    } else if( !n1->m_data->isTip && n2->m_data->isTip ) {
        cont.push_back( rooted_bifurcation<lnode>( n, n2, n1, rooted_bifurcation<lnode>::TIP_INNER ));
        rooted_traveral_order_rec( n1, cont );
    } else {
        cont.push_back( rooted_bifurcation<lnode>( n, n1, n2, rooted_bifurcation<lnode>::INNER_INNER ));
        rooted_traveral_order_rec( n1, cont );
        rooted_traveral_order_rec( n2, cont );
    }
}

template <class lnode, class container>
void rooted_traveral_order( lnode *n1, lnode *n2, container &cont ) {
    
    if( !n1->m_data->isTip ) {
        rooted_traveral_order_rec<lnode, container>( n1, cont );
    }
    if( !n2->m_data->isTip ) {
        rooted_traveral_order_rec<lnode, container>( n2, cont );
    }
    
    
    std::reverse( cont.begin(), cont.end());
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

#endif