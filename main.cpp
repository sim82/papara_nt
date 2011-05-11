#include "ivymike/multiple_alignment.h"
#include <stdexcept>
#include <iostream>
#include <vector>
#include <deque>
#include <map>

#include "parsimony.h"
#include "pars_align_vec.h"
#include "fasta.h"


#include "ivymike/tree_parser.h"
#include "ivymike/time.h"


using namespace ivy_mike;
using namespace ivy_mike::tree_parser_ms;

class ostream_test {
    std::ostream &m_os;
    
public:
    ostream_test( std::ostream &os ) : m_os(os) {}
    void operator()(int i) {
        m_os << i;
    }
};

typedef unsigned int parsimony_state;

static inline parsimony_state dna_to_parsimony_state( uint8_t c ) {
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
        
        case 'T':
        case 't':
            return 0x8;
            
        default:
            return 0xf;
    };
    
}


static inline int dna_to_cgap( uint8_t c ) {
    switch( c ) {
    case 'A':
    case 'a':
    case 'C':
    case 'c':
    case 'G':
    case 'g':
    case 'T':
    case 't':
        return 0x0;
        
    default:
        return AUX_CGAP;
    };
}

class pvec_cgap {
//     aligned_buffer<parsimony_state> v;
      std::vector<parsimony_state> v;
      std::vector<int> auxv;
      
public:
    void init( const std::vector<uint8_t> &seq ) {
        assert( v.size() == 0 );
        v.resize(seq.size());
        auxv.resize( seq.size() );
        std::transform( seq.begin(), seq.end(), v.begin(), dna_to_parsimony_state );
        std::transform( seq.begin(), seq.end(), auxv.begin(), dna_to_cgap );
    }
    
    static void newview( pvec_cgap &p, pvec_cgap &c1, pvec_cgap &c2, tip_case tc ) {
        assert( c1.v.size() == c2.v.size() );
        
//         p.v.resize(0);
        p.v.resize(c1.v.size());
        p.auxv.resize(c1.auxv.size());
        
        
        for( size_t i = 0; i < c1.v.size(); i++ ) {
            parsimony_state ps = c1.v[i] & c2.v[i];
            
            if( ps == 0 ) {
                ps = c1.v[i] | c2.v[i];
            }
            
            //p.v.push_back( ps );
            p.v[i] = ps;
            
        
            
            const int a1 = c1.auxv[i];
            const int a2 = c2.auxv[i];
            
            const bool cgap1 = (a1 & AUX_CGAP) != 0;
            const bool cgap2 = (a2 & AUX_CGAP) != 0;
            
//             const bool open1 = (a1 & AUX_OPEN) != 0;
//             const bool open2 = (a2 & AUX_OPEN) != 0;
            
            p.auxv[i] = 0;
            
            if( tc == TIP_TIP ) {
                if( cgap1 && cgap2 ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( cgap1 != cgap2 ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            } else if( tc == TIP_INNER ) {
                if( cgap1 && cgap2 ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( cgap1 != cgap2 ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            } else {
                if( a1 == AUX_CGAP && a2 == AUX_CGAP ) {
                    p.auxv[i] = AUX_CGAP;
                } else if( a1 == AUX_CGAP || a2 == AUX_CGAP ) {
                    p.auxv[i] = AUX_CGAP | AUX_OPEN;
                }
            }

        }
    }
    
    inline size_t size() {
        return v.size();
    }
    
    inline void to_int_vec( std::vector<int> &outv ) {
        outv.resize( v.size() );
        
        std::copy( v.begin(), v.end(), outv.begin() );
    }
    
    inline void to_aux_vec( std::vector<unsigned int> &outv ) {
        outv.resize( v.size() );
        std::copy( auxv.begin(), auxv.end(), outv.begin() );
        
         
//         std::for_each( auxv.begin(), auxv.end(), ostream_test(std::cout) );
        
        
    }
  
  
};


// class pvec_ugly {
//     parsimony_state *m_v;
//     size_t m_size;
//     
// public:
//     pvec_ugly() : m_v(0), m_size(0) {}
//     ~pvec_ugly() {
//         delete[] m_v;
//     }
//     void init( const std::vector<uint8_t> &seq ) {
//         assert( m_v == 0 );
//         delete[] m_v;
//         m_size = seq.size(0);
//         m_v = new parsimony_state[m_size];
//         
//         std::transform( seq.begin(), seq.end(), m_v, dna_to_parsimony_state );
//     }
//     
//     static void newview( pvec_ugly &p, pvec_ugly &c1, pvec_ugly &c2 ) {
//         assert( c1.m_size() == c2.m_size() );
//         
// //         p.v.resize(0);
//         //p.v.resize(c1.v.size());
// #error continue here!
//         
//         for( size_t i = 0; i < c1.v.size(); i++ ) {
//             parsimony_state ps = c1.v[i] & c2.v[i];
//             
//             if( ps == 0 ) {
//                 ps = c1.v[i] | c2.v[i];
//             }
//             
//             //p.v.push_back( ps );
//             p.v[i] = ps;
//         
//         }
//     }
//     
// };


template<class pvec_t>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {
//     static int ct;
    //std::vector<parsimony_state> m_pvec;
    pvec_t m_pvec;
    
public:
//     int m_ct;
    my_adata_gen() {
     
//         std::cout << "my_adata\n";
        
    }
        
    virtual ~my_adata_gen() {
     
//         std::cout << "~my_adata\n";
        
    }
    
    virtual void visit() {
//         std::cout << "tr: " << m_ct << "\n";
    }
    void init_pvec(const std::vector< uint8_t >& seq) {
        
        
        m_pvec.init( seq );
//         std::cout << "init_pvec: " << m_pvec.size() << "\n";
//                 m_pvec.reserve(seq.size());
//         for( std::vector< uint8_t >::const_iterator it = seq.begin(); it != seq.end(); ++it ) {
//             m_pvec.push_back(dna_to_parsimony_state(*it));
//             
//         }
    }
    pvec_t &get_pvec() {
        return m_pvec;
    }
    
 
};

// inline void newview_parsimony( std::vector<parsimony_state> &p, const std::vector<parsimony_state> &c1, const std::vector<parsimony_state> &c2 ) {
//     
// }



// inline std::ostream &operator<<( std::ostream &os, const my_adata &rb ) {
// 
//     os << "my_adata: " << rb.m_ct;
// }

template<class ndata_t>
class my_fact : public ivy_mike::tree_parser_ms::node_data_factory {
  
    virtual ndata_t *alloc_adata() {
     
        return new ndata_t;
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

// template<class lnode>
// void traverse( lnode *n ) {
//     n->m_data->visit();
//
//     if( n->back != 0 ) {
//         traverse_rec(n->back);
//     }
//
//     if( n->next->back != 0 ) {
//         traverse_rec(n->next->back);
//     }
//
//     if( n->next->next->back != 0 ) {
//         traverse_rec(n->next->next->back);
//     }
//
//
// }

template<class pvec_t>
void do_newview( pvec_t &root_pvec, lnode *n1, lnode *n2, bool incremental ) {
    typedef my_adata_gen<pvec_t> my_adata;
    
    std::deque<rooted_bifurcation<lnode> > trav_order;
    
    std::cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";
    
    rooted_traveral_order( n1, n2, trav_order, incremental );
    
    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
//         std::cout << *it << "\n";
        
        my_adata *p = dynamic_cast<my_adata *>( it->parent->m_data.get());
        my_adata *c1 = dynamic_cast<my_adata *>( it->child1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( it->child2->m_data.get());
//         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;
        
        
        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->tc);
        
    }
    
    
    
    
    
    {
        my_adata *c1 = dynamic_cast<my_adata *>( n1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( n2->m_data.get());
        
//         tip_case tc;
        
        if( c1->isTip && c2->isTip ) {
        
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), TIP_TIP );
        } else if( c1->isTip && !c2->isTip ) {
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), TIP_INNER );
        } else if( !c1->isTip && c2->isTip ) {
            pvec_t::newview(root_pvec, c2->get_pvec(), c1->get_pvec(), TIP_INNER );
        } else {
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), INNER_INNER );
        }
        
        
    }
//     std::cout << std::hex;
//     for( std::vector< parsimony_state >::const_iterator it = root_pvec.begin(); it != root_pvec.end(); ++it ) {
//         std::cout << *it;
//     }
//     
//     std::cout << std::dec << std::endl;
    
}


static void seq_to_nongappy_pvec( std::vector<uint8_t> &seq, std::vector<uint8_t> &pvec ) {
    pvec.resize( 0 );
    
    for( unsigned int i = 0; i < seq.size(); i++ ) {
        uint8_t ps;
        
        switch( seq[i] ) {
        case 'A':
        case 'a':
            ps = 0x1;
            break;
        case 'C':
        case 'c':
            ps = 0x2;
            break;
            
        case 'G':
        case 'g':
            ps = 0x4;
            break;
            
        case 'T':
        case 't':
            ps = 0x8;
            break;
            
        default:
            ps = 0;
        }
        
        if( ps != 0 ) {
            pvec.push_back(ps);
        }
        
    }
    
}

void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq );


int mainx() {

    
    
    typedef pvec_cgap pvec_t;
	
    typedef my_adata_gen<pvec_t> my_adata;
    
    ivy_mike::timer t;
    ivy_mike::tree_parser_ms::ln_pool pool( sptr::shared_ptr<my_fact<my_adata> >( new my_fact<my_adata> ) );
    ivy_mike::tree_parser_ms::parser tp( "test_1604/RAxML_bestTree.ref_orig", pool );
    ivy_mike::tree_parser_ms::lnode * n = tp.parse();
    
    n = towards_tree( n );
    
    
    
    
    
    typedef tip_collector<lnode> tc_t;
    
    tc_t tc;
    
    visit_lnode( n, tc );
    
    std::map<std::string, sptr::shared_ptr<lnode> > name_to_lnode;
    
    for( std::vector< sptr::shared_ptr<lnode> >::iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
        // 		std::cout << (*it)->m_data->tipName << "\n";
        name_to_lnode[(*it)->m_data->tipName] = *it;
    }
    
    multiple_alignment ma;
    ma.load_phylip( "test_1604/orig.phy.1" );
    for( unsigned int i = 0; i < ma.names.size(); i++ ) {
        
        
        sptr::shared_ptr< lnode > ln = name_to_lnode[ma.names[i]];
        // 		adata *ad = ln->m_data.get();
        
        assert( typeid(*ln->m_data.get()) == typeid(my_adata ) );
        my_adata *adata = static_cast<my_adata *> (ln->m_data.get());
        
        adata->init_pvec( ma.data[i] );
        
    }
    
    
    //     traverse( n );
    
    edge_collector<lnode> ec;
    visit_edges( n, ec );
    
    std::cout << "edges: " << ec.m_edges.size() << "\n";
    
    std::vector< pvec_t > m_parsvecs;
    m_parsvecs.resize( ec.m_edges.size() );
    
    
//    mapped_file qsf( "test_1604/qs.fa" );
	std::ifstream qsf( "test_1604/qs.fa" );
    std::vector<std::string> qs_names;
    std::vector<std::vector<uint8_t> > qs_seqs;
    
    std::vector<std::vector <uint8_t> > qs_nongappy;
    
    
    
    
    read_fasta( qsf, qs_names, qs_seqs);
    
    
    
    qs_nongappy.resize( qs_names.size() );
    std::vector <int> qs_bestscore(qs_names.size());
    std::fill( qs_bestscore.begin(), qs_bestscore.end(), 32000);
    std::vector <int> qs_bestedge(qs_names.size());
    
    const size_t VW = pars_align_vec::WIDTH;
    pars_align_vec::arrays<VW> arrays;
    
    int n_groups = (ec.m_edges.size() / VW) + 1;
    
    std::vector<int> seqlist[VW];
    const int *seqptrs[VW];
    std::vector<unsigned int> auxlist[VW];
    const unsigned int *auxptrs[VW];
    pvec_t root_pvec;
    
    for ( int j = 0; j < n_groups; j++ ) {
        int num_valid = 0;
    
        for( unsigned int i = 0; i < VW; i++ ) {
            
            unsigned int edge = j * VW + i;
            if( edge < ec.m_edges.size() ) {
                
                
                do_newview( root_pvec, ec.m_edges[edge].first, ec.m_edges[edge].second, true );
                root_pvec.to_int_vec(seqlist[i]);
                root_pvec.to_aux_vec(auxlist[i]);
                
                seqptrs[i] = seqlist[i].data();
                auxptrs[i] = auxlist[i].data();
                num_valid++;
            } else {
                if( i < 1 ) {
                    throw std::runtime_error( "bad integer mathematics" );
                }
                seqlist[i] = seqlist[i-1];
            }
            
        }
        
        
        
        for( unsigned int i = 0; i < qs_names.size(); i++ ) {
            if( qs_nongappy[i].size() == 0 ) {
                seq_to_nongappy_pvec( qs_seqs[i], qs_nongappy[i] );    
            }
            size_t stride = 1;
            size_t aux_stride = 1;
            pars_align_vec pa( seqptrs, qs_nongappy[i].data(), seqlist[0].size(), qs_nongappy[i].size(), stride, auxptrs, aux_stride, arrays, 0 );
            
            
            pars_align_vec::score_t *score_vec = pa.align_freeshift();
            for( int k = 0; k < num_valid; k++ ) {
                if( score_vec[k] < qs_bestscore[i] ) {
                    qs_bestscore[i] = score_vec[k];
                    qs_bestedge[i] = j * VW + k;
                }
            }
        }
     
    }
    
    for( unsigned int i = 0; i < qs_names.size(); i++ ) {
        std::cout << qs_names[i] << " " << qs_bestedge[i] << " " << qs_bestscore[i] << "\n";
        
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


    std::cout << t.elapsed() << std::endl;
    return 0;
//     getchar();
}


int main() {
    mainx();
    return 0;
}