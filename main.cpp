#include <ivymike/MultipleAlignment.h>

#include <iostream>
#include <vector>
#include <deque>
#include <map>

#include "parsimony.h"
#include "ivymike/tree_parser.h"
#include "ivymike/time.h"

typedef unsigned int parsimony_state;

class my_adata : public ivy_mike::tree_parser_ms::adata {
    static int ct;
    std::vector<parsimony_state> m_pvec;
    
public:
    int m_ct;
    my_adata() : m_ct(ct++) {
     
//         std::cout << "my_adata\n";
        
    }
        
    virtual ~my_adata() {
     
//         std::cout << "~my_adata\n";
        
    }
    
    virtual void visit() {
        std::cout << "tr: " << m_ct << "\n";
    }
};

int my_adata::ct = 0;

// inline std::ostream &operator<<( std::ostream &os, const my_adata &rb ) {
// 
//     os << "my_adata: " << rb.m_ct;
// }


class my_fact : public ivy_mike::tree_parser_ms::node_data_factory {
  
    virtual my_adata *alloc_adata() {
     
        return new my_adata;
    }
    
};



template<class lnode>
void traverse_rec( lnode *n ) {

    n->m_data->visit();
    
    if( n->next->back != 0 ) {
        traverse_rec(n->next->back);    
    }
    
    if( n->next->next->back != 0 ) {
        traverse_rec(n->next->next->back);    
    }
}

template<class lnode>
void traverse( lnode *n ) {
    n->m_data->visit();
    
    if( n->back != 0 ) {
        traverse_rec(n->back);    
    }
    
    if( n->next->back != 0 ) {
        traverse_rec(n->next->back);    
    }
    
    if( n->next->next->back != 0 ) {
        traverse_rec(n->next->next->back);    
    }
    
    
}





int main() {
//     getchar();
    //ivymike::TreeParser tp( "./RAxML_bipartitions.1604.BEST.WITH" );
    using namespace ivy_mike::tree_parser_ms;
	
    ivy_mike::timer t;
    ivy_mike::tree_parser_ms::ln_pool pool( boost::shared_ptr<my_fact>( new my_fact ) );
    ivy_mike::tree_parser_ms::parser tp( "small.tree", pool );
    ivy_mike::tree_parser_ms::lnode * n = tp.parse();
    
    n = towards_tree( n );
    
    
    
    traverse<ivy_mike::tree_parser_ms::lnode>( n );
    
    std::deque<rooted_bifurcation<ivy_mike::tree_parser_ms::lnode> > trav_order;
    
    std::cout << "traversal for branch: " << *(n->m_data) << " " << *(n->back->m_data) << "\n";
    
    rooted_traveral_order( n, n->back, trav_order );
    
    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
        std::cout << *it << "\n";
    }
    
    typedef tip_collector<lnode> tc_t;
	
	tc_t tc;
	
	visit_lnode< tip_collector<ivy_mike::tree_parser_ms::lnode> >( n, tc );

	std::map<std::string, boost::shared_ptr<lnode> > name_to_lnode;
	
	for( std::vector< ivy_mike::tree_parser_ms::lnode* >::iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
		std::cout << (*it)->m_data->tipName << "\n";
		name_to_lnode[(*it)->m_data->tipName] = boost::shared_ptr<lnode>((*it)->get_smart_ptr());
	}
    
    printf( "n: %f %d\n", n->backLen, n->m_data->isTip );
    
    assert( n->next->back != 0 && n->next->next->back != 0 );
    
    n->next->back->back = n->next->next->back;
    n->next->next->back->back = n->next->back;
    n = n->next->back;
    
    
    {
        ivy_mike::timer t2;
        pool.clear();
        pool.mark(n);
        pool.sweep();
        
        std::cout << t2.elapsed() << std::endl;
    }
    
    
    {
        ivy_mike::timer t2;
        pool.clear();
     //   pool.mark(n);
        pool.sweep();
        
        std::cout << t2.elapsed() << std::endl;
    }
    //ivymike::LN *n = tp.parse();
    
//     getchar();
    //ivymike::LN::free( n );
//     delete n;
//     getchar();
    
    std::cout << t.elapsed() << std::endl;
}