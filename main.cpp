#include "ivymike/multiple_alignment.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <map>
#include <cstring>
#include "parsimony.h"
#include "pars_align_vec.h"
#include "pars_align_seq.h"
#include "fasta.h"


#include "ivymike/tree_parser.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/thread.h"


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

typedef int parsimony_state;


// this class is so strange, we could name a new design-pattern after it.
// it does something like a static {} block in java...
class dna_parsimony_mapping_real {
    std::vector<uint8_t> m_p2d;
    std::vector<parsimony_state> m_d2p;

public:
    dna_parsimony_mapping_real() : m_p2d(16), m_d2p(256, -1) 
    {
        const uint8_t pd[16] = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'};
        m_p2d.assign( pd, pd + 16 );
        
        for( int i = 0; i < 16; i++ ) {
            m_d2p[pd[i]] = i;
            m_d2p[std::tolower(pd[i])] = i;
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
            throw std::runtime_error( "illegal dna character" );
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

// static inline parsimony_state dna_to_parsimony_state( uint8_t c ) {
//     switch( c ) {
//         case 'A':
//         case 'a':
//             return 0x1;
//             
//         case 'C':
//         case 'c':
//             return 0x2;
//             
//         case 'G':
//         case 'g':
//             return 0x4;
//         
//         case 'U':
//         case 'u':
//         case 'T':
//         case 't':
//             return 0x8;
//             
//         default:
//             return 0xf;
//     };
//     
// }
// 
// 
// static inline int dna_to_cgap( uint8_t c ) {
//     switch( c ) {
//     case 'A':
//     case 'a':
//     case 'C':
//     case 'c':
//     case 'G':
//     case 'g':
//     case 'T':
//     case 't':
//     case 'U':
//     case 'u':
//         return 0x0;
//         
//     default:
//         return AUX_CGAP;
//     };
// }

static bool g_dump_aux = false;

class pvec_cgap {
    //     aligned_buffer<parsimony_state> v;
    std::vector<parsimony_state> v;
    std::vector<int> auxv;
    
public:
    void init( const std::vector<uint8_t> &seq ) {
        assert( v.size() == 0 );
        v.resize(seq.size());
        auxv.resize( seq.size() );
        std::transform( seq.begin(), seq.end(), v.begin(), dna_parsimony_mapping::d2p );
        std::transform( seq.begin(), seq.end(), auxv.begin(), dna_parsimony_mapping::d2aux );
    }
    
    static void newview( pvec_cgap &p, pvec_cgap &c1, pvec_cgap &c2, tip_case tc ) {
        assert( c1.v.size() == c2.v.size() );
        
//         p.v.resize(0);
        p.v.resize(c1.v.size());
        p.auxv.resize(c1.auxv.size());
        
        if( g_dump_aux ) {
            std::cout << "1:";
            std::copy( c1.auxv.begin(), c1.auxv.end(), std::ostream_iterator<int>(std::cout) );
            std::cout << "\n2:";
            std::copy( c2.auxv.begin(), c2.auxv.end(), std::ostream_iterator<int>(std::cout) );
            std::cout << "\n";
        }
        
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
//         std::cout << "v: " << v.size() << "\n";

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
    
    //std::cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";
    
    rooted_traveral_order( n1, n2, trav_order, incremental );
//     std::cout << "traversal: " << trav_order.size() << "\n";
    
    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
//         std::cout << *it << "\n";
        
        my_adata *p = dynamic_cast<my_adata *>( it->parent->m_data.get());
        my_adata *c1 = dynamic_cast<my_adata *>( it->child1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( it->child2->m_data.get());
//         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;
        
//         std::cout << "tip case: " << (*it) << "\n";
        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->tc);
        
    }
    
    
    
    
    
    {
        my_adata *c1 = dynamic_cast<my_adata *>( n1->m_data.get());
        my_adata *c2 = dynamic_cast<my_adata *>( n2->m_data.get());
        
//         tip_case tc;
        
        if( c1->isTip && c2->isTip ) {
//                 std::cout << "root: TIP TIP\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), TIP_TIP );
        } else if( c1->isTip && !c2->isTip ) {
//                 std::cout << "root: TIP INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), TIP_INNER );
        } else if( !c1->isTip && c2->isTip ) {
//                 std::cout << "root: INNER TIP\n";
            pvec_t::newview(root_pvec, c2->get_pvec(), c1->get_pvec(), TIP_INNER );
        } else {
//                 std::cout << "root: INNER INNER\n";
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
        uint8_t ps = dna_parsimony_mapping::d2p(seq[i]);
  
        if( ps == 0x1 || ps == 0x2 || ps == 0x4 || ps == 0x8 ) {
            pvec.push_back(ps);
        }
        
    }
    
}

void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq );




class papara_nt {
    
    typedef pvec_cgap pvec_t;
    typedef my_adata_gen<pvec_t> my_adata;
    
        
    const static size_t VW = pars_align_vec::WIDTH;
    
    struct block_t {
        block_t() {
            memset( this, 0, sizeof( block_t )); // FIXME: hmm, this is still legal?
        }
        
        // WARNING: these are pointers into m_ref_pvecs and m_ref_aux
        // make sure they stay valid!
        const int *seqptrs[VW];
        const unsigned int *auxptrs[VW];
        size_t ref_len;
        int edges[VW];
        int num_valid;
    };
    
    papara_nt( const papara_nt &other );
    papara_nt & operator=( const papara_nt &other );
    
    
    
    multiple_alignment m_ref_ma;
    std::auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    edge_collector<lnode> m_ec;
    
    std::vector<std::string> m_qs_names;
    std::vector<std::vector<uint8_t> > m_qs_seqs;
        
    std::vector<std::vector <uint8_t> > m_qs_pvecs;
    
    std::vector<std::vector <int> > m_ref_pvecs;
    std::vector<std::vector <unsigned int> > m_ref_aux;
    
    ivy_mike::mutex m_qmtx; // mutex for the block queue and the qs best score/edge arrays
    std::deque<block_t> m_blockqueue;
    std::vector <int> m_qs_bestscore;
    std::vector <int> m_qs_bestedge;
        
    
    
    
    class worker {
        papara_nt &m_pnt;
    public:
        worker( papara_nt & pnt ) : m_pnt(pnt) {}
        void operator()() {
            
            pars_align_vec::arrays<VW> arrays;
            pars_align_seq::arrays seq_arrays(true);
            
            
            while( true ) {
                block_t block;
                
                {
                    ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );
                    if( m_pnt.m_blockqueue.empty() ) {
                        break;
                    }
                    block = m_pnt.m_blockqueue.front();
                    m_pnt.m_blockqueue.pop_front();
                }
                
                
                
                
                for( unsigned int i = 0; i < m_pnt.m_qs_names.size(); i++ ) {
                        
                    size_t stride = 1;
                    size_t aux_stride = 1;
                    pars_align_vec pa( block.seqptrs, m_pnt.m_qs_pvecs[i].data(), block.ref_len, m_pnt.m_qs_pvecs[i].size(), stride, block.auxptrs, aux_stride, arrays, 0 );
                    
                    
                    pars_align_vec::score_t *score_vec = pa.align_freeshift();
                
                    {
                        ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );
                        
                        for( int k = 0; k < block.num_valid; k++ ) {
                            
                            
                            
                            if( score_vec[k] < m_pnt.m_qs_bestscore[i] || (score_vec[k] == m_pnt.m_qs_bestscore[i] && block.edges[k] < m_pnt.m_qs_bestedge[i] )) {
                                const bool validate = false;
                                if( validate ) {
                                    const int *seqptr = block.seqptrs[k];
                                    const unsigned int *auxptr = block.auxptrs[k];
                                    
                                    pars_align_seq pas( seqptr, m_pnt.m_qs_pvecs[i].data(), block.ref_len, m_pnt.m_qs_pvecs[i].size(), stride, auxptr, aux_stride, seq_arrays, 0 );
                                    int res = pas.alignFreeshift(INT_MAX);
                                    
                                    if( res != score_vec[k] ) {
                                        
                                    
                                        std::cout << "meeeeeeep! score: " << score_vec[k] << " " << res << "\n";
                                    }
                                }
                                
                                m_pnt.m_qs_bestscore[i] = score_vec[k];
                                m_pnt.m_qs_bestedge[i] = block.edges[k];
                            }
                        }
                    }
                }
                
            }
            
            
        }
    };
    
    void build_block_queue() {
        // creates the list of ref-block to be consumed by the worker threads.  A ref-block onsists of N ancestral state sequences, where N='width of the vector unit'.
        // The vectorized alignment implementation will align a QS against a whole ref-block at a time, rather than a single ancestral state sequence as in the 
        // sequencial algorithm. 
        
        assert( m_blockqueue.empty() );
        int n_groups = (m_ec.m_edges.size() / VW);
        if( (m_ec.m_edges.size() % VW) != 0 ) {
            n_groups++;
        }
        
        
//         std::vector<int> seqlist[VW];
//         const int *seqptrs[VW];
//         std::vector<unsigned int> auxlist[VW];
//         const unsigned int *auxptrs[VW];
        
        
        
        for ( int j = 0; j < n_groups; j++ ) {
            int num_valid = 0;
            
            size_t ref_size = 0;
            
            block_t block;
            
            for( unsigned int i = 0; i < VW; i++ ) {
            
                unsigned int edge = j * VW + i;
                if( edge < m_ec.m_edges.size()) {
                    block.edges[i] = edge;
                    block.num_valid++;
                    
                    block.seqptrs[i] = m_ref_pvecs[edge].data();
                    block.auxptrs[i] = m_ref_aux[edge].data();
                    block.ref_len = m_ref_pvecs[edge].size();
                    //                     do_newview( root_pvec, m_ec.m_edges[edge].first, m_ec.m_edges[edge].second, true );
//                     root_pvec.to_int_vec(seqlist[i]);
//                     root_pvec.to_aux_vec(auxlist[i]);
//                     
//                     seqptrs[i] = seqlist[i].data();
//                     auxptrs[i] = auxlist[i].data();\

                    num_valid++;
                } else {
                    if( i < 1 ) {
                        std::cout << "edge: " << edge << " " << m_ec.m_edges.size() << std::endl;
                        
                        throw std::runtime_error( "bad integer mathematics" );
                    }
                    block.edges[i] = block.edges[i-1];
                    
                    block.seqptrs[i] = block.seqptrs[i-1];
                    block.auxptrs[i] = block.auxptrs[i-1];
                }
                
            }
            m_blockqueue.push_back(block);
        }   
    }
    
    void build_ref_vecs() {
        // pre-create the ancestral state vectors. This step is necessary for the threaded version, because otherwise, each
        // thread would need an independent copy of the tree to do concurrent newviews. Anyway, having a copy of the tree
        // in each thread will most likely use more memory than storing the pre-calculated vectors. 
        
        // TODO: maybe try lazy create/cache of the asv's in the threads
        
        assert( m_ref_aux.empty() && m_ref_pvecs.empty() );
        
        m_ref_pvecs.reserve( m_ec.m_edges.size() );
        m_ref_aux.reserve( m_ec.m_edges.size() );
        
        
        for( size_t i = 0; i < m_ec.m_edges.size(); i++ ) {
            pvec_t root_pvec;
            
            std::cout << "newview for branch " << i << ": " << *(m_ec.m_edges[i].first->m_data) << " " << *(m_ec.m_edges[i].second->m_data) << "\n";
            
            if( i == 340 ) {
                g_dump_aux = true;
            }
            
            do_newview( root_pvec, m_ec.m_edges[i].first, m_ec.m_edges[i].second, true );
            
            g_dump_aux = false;
            // TODO: try something fancy with rvalue refs...
            
            m_ref_pvecs.push_back( std::vector<int>() );
            m_ref_aux.push_back( std::vector<unsigned int>() );
            
            root_pvec.to_int_vec(m_ref_pvecs.back());
            root_pvec.to_aux_vec(m_ref_aux.back());
                    
            
        }
           
    }
    
public:    
    
    
    
    papara_nt( const char* opt_tree_name, const char *opt_alignment_name, const char *opt_qs_name ) 
      : m_ln_pool(new ivy_mike::tree_parser_ms::ln_pool( sptr::shared_ptr<my_fact<my_adata> >( new my_fact<my_adata> ))) 
    {
     
        // load input data: ref-tree, ref-alignment and query sequences
        
        //
        // parse the reference tree
        //
        
        
        ln_pool &pool = *m_ln_pool;
        tree_parser_ms::parser tp( opt_tree_name, pool );
        tree_parser_ms::lnode * n = tp.parse();
        
        n = towards_tree( n );
        //
        // create map from tip names to tip nodes
        //
        typedef tip_collector<lnode> tc_t;
        tc_t tc;
        
        visit_lnode( n, tc );
        
        std::map<std::string, sptr::shared_ptr<lnode> > name_to_lnode;
        
        for( std::vector< sptr::shared_ptr<lnode> >::iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
            std::cout << (*it)->m_data->tipName << "\n";
            name_to_lnode[(*it)->m_data->tipName] = *it;
        }
        
        //
        // read reference alignment: store the ref-seqs in the tips of the ref-tree
        //
        m_ref_ma.load_phylip( opt_alignment_name );
        for( unsigned int i = 0; i < m_ref_ma.names.size(); i++ ) {
            
            
            sptr::shared_ptr< lnode > ln = name_to_lnode[m_ref_ma.names[i]];
            //      adata *ad = ln->m_data.get();
            
            assert( typeid(*ln->m_data.get()) == typeid(my_adata ) );
            my_adata *adata = static_cast<my_adata *> (ln->m_data.get());
            
            adata->init_pvec( m_ref_ma.data[i] );
            
        }
        
        //
        // collect list of edges
        //
        
        visit_edges( n, m_ec );
        
        std::cout << "edges: " << m_ec.m_edges.size() << "\n";
        
//         std::vector< pvec_t > m_parsvecs;
//         m_parsvecs.resize( m_ec.m_edges.size() );
        
        
        //
        // read query sequences
        //
        std::ifstream qsf( opt_qs_name );
        
        
        read_fasta( qsf, m_qs_names, m_qs_seqs);
        
        //
        // setup qs best-score/best-edge lists
        //
        
        
        m_qs_pvecs.resize( m_qs_names.size() );
        
        m_qs_bestscore.resize(m_qs_names.size());
        std::fill( m_qs_bestscore.begin(), m_qs_bestscore.end(), 32000);
        m_qs_bestedge.resize(m_qs_names.size());
        
        //
        // pre-create the reference pvecs/auxvecs
        //
        
        build_ref_vecs();
          
        
        //
        // preprocess query sequences
        //
        
        for( int i = 0; i < m_qs_seqs.size(); i++ ) {
            seq_to_nongappy_pvec( m_qs_seqs[i], m_qs_pvecs[i] );    
        }
        
     
    }
    
    
    void calc_scores( const size_t n_threads ) {
        
        //
        // build the alignment blocks
        //
           
        
        build_block_queue();
        
        //
        // work
        //
        ivy_mike::timer t1;
        ivy_mike::thread_group tg;
        
        while( tg.size() < n_threads ) {
            tg.create_thread(worker(*this));
        }
        tg.join_all();
        
        std::cerr << "scoring finished: " << t1.elapsed() << "\n";
           
    }
    
    void print_best_scores( std::ostream &os ) {
        
        os << std::setfill ('0');
        for( unsigned int i = 0; i < m_qs_names.size(); i++ ) {
            os << m_qs_names[i] << " "  << std::setw (4) << m_qs_bestedge[i] << " " << std::setw(5) << m_qs_bestscore[i] << "\n";
            
        }   
    }
    
    
    void gapstream_to_alignment( const std::vector<uint8_t> &gaps, const std::vector<uint8_t> &raw, std::vector<uint8_t> &out, uint8_t gap_char ) {
        
        std::vector<uint8_t>::const_reverse_iterator rit = raw.rbegin();

        for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {

            if ( *git == 1) {
                out.push_back(gap_char);
            } else if ( *git == 0 ) {

                out.push_back(*rit);
                ++rit;
            } else {
                ++rit; // just consume one QS character
            }
        }

        std::reverse(out.begin(), out.end());
    }
    
    
    static uint8_t normalize_dna( uint8_t c ) {
        c = std::toupper(c);
        
        if( c == 'U' ) {
            c = 'T';
        }
        
        return c;
    }
    
    static char num_to_ascii( int n ) {
        if( n >= 0 && n <= 9 ) {
            return '0' + n;
        } else if( n >= 0xa && n <= 0xf ) {
            return 'a' + n;
        } else {
            throw std::runtime_error( "not a single digit (hex) number" );
        }
    }
    
    void align_best_scores( std::ostream &os ) {
        // create the actual alignments for the best scoring insertion position (=do the traceback)
        
        pars_align_seq::arrays seq_arrays(true);
        
        for( unsigned int i = 0; i < m_qs_names.size(); i++ ) {
            int best_edge = m_qs_bestedge[i];
            
            assert( best_edge >= 0 && best_edge < m_ref_pvecs.size() );
            
            const int *seqptr = m_ref_pvecs[best_edge].data();
            const unsigned int *auxptr = m_ref_aux[best_edge].data();
            
            const size_t ref_len = m_ref_pvecs[best_edge].size();
            
            const size_t stride = 1;
            const size_t aux_stride = 1;
            pars_align_seq pas( seqptr, m_qs_pvecs[i].data(), ref_len, m_qs_pvecs[i].size(), stride, auxptr, aux_stride, seq_arrays, 0 );
            int res = pas.alignFreeshift(INT_MAX);
            
            
            std::vector<uint8_t> tbv;
            pas.tracebackCompressed(tbv);
            
            std::vector<uint8_t> out_qs;
            
            
            gapstream_to_alignment(tbv, m_qs_pvecs[i], out_qs, 0xf);
            
            std::string out_string;
            out_string.resize(out_qs.size() );
            
            std::transform( out_qs.begin(), out_qs.end(), out_string.begin(), dna_parsimony_mapping::p2d );
            
            os << m_qs_names[i] << "\t" << out_string << "\n";
            
            const bool dump_auxvec = false;
            if( dump_auxvec )
            {
                std::string auxv;
                auxv.resize(m_ref_aux[best_edge].size());
                
                std::transform( m_ref_aux[best_edge].begin(), m_ref_aux[best_edge].end(), auxv.begin(), num_to_ascii );
                os << m_qs_names[i] << "\t" << auxv << "\n";
            }
            
            
            if( res != m_qs_bestscore[i] ) {
                std::cout << "meeeeeeep! score: " << m_qs_bestscore[i] << " " << res << "\n";
            }
        }
    }
    
    ~papara_nt() {
        
        ivy_mike::timer t2;
        m_ln_pool->clear();
        //   pool.mark(n);
        m_ln_pool->sweep();
        
        std::cout << t2.elapsed() << std::endl;

    }
    void dump_ref_seqs ( std::ostream &os ) {
        for( size_t i = 0; i < m_ref_ma.names.size(); i++ ) {
            std::string outs;
            outs.resize(m_ref_ma.data[i].size() );
            
            std::transform( m_ref_ma.data[i].begin(), m_ref_ma.data[i].end(), outs.begin(), normalize_dna );
            
            os << m_ref_ma.names[i] << "\t" << outs << "\n";
        }
    }
    
    
    void write_result_phylip( std::ostream &os ) {
        os << " 1 2\n";
        dump_ref_seqs(os);
        align_best_scores(os);
        
    }
};

int main( int argc, char *argv[] ) {

    namespace igo = ivy_mike::getopt;
    
    ivy_mike::getopt::parser igp;
    
    std::string opt_tree_name;
    std::string opt_alignment_name;
    std::string opt_qs_name;
    igp.add_opt( 't', igo::value<std::string>(opt_tree_name) );
    igp.add_opt( 's', igo::value<std::string>(opt_alignment_name) );
    igp.add_opt( 'q', igo::value<std::string>(opt_qs_name) );
    
    igp.parse(argc,argv);
    
    if( igp.opt_count('t') != 1 || igp.opt_count('s') != 1 || igp.opt_count('q') != 1 ) {
        std::cerr << "missing options -t, -q and/or -s\n";
        return 0;
    }
    ivy_mike::timer t;
    
    papara_nt pnt( opt_tree_name.c_str(), opt_alignment_name.c_str(), opt_qs_name.c_str() );
    pnt.calc_scores( 4 );
    
    {
        std::ofstream os( "papara_scores.txt" );
        pnt.print_best_scores(os);
    }
    
    {
        std::ofstream os( "papara_ali.phy" );
//         pnt.dump_ref_seqs(os);
//         pnt.align_best_scores(os);
        pnt.write_result_phylip(os);
    }
    
    
    
    

    //ivymike::LN *n = tp.parse();
    
//     getchar();
    //ivymike::LN::free( n );
//     delete n;


    std::cout << t.elapsed() << std::endl;
    return 0;
//     getchar();
}

