
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <map>
#include <functional>
#include <cstring>
#include <boost/io/ios_state.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>



#include "sequence_model.h"
#include "parsimony.h"
#include "pvec.h"
#include "align_utils.h"

#include "pars_align_seq.h"
#include "pars_align_gapp_seq.h"
#include "fasta.h"
#include "vec_unit.h"
#include "align_pvec_vec.h"
#include "align_utils.h"

#include "ivymike/tree_parser.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/thread.h"
#include "ivymike/demangle.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
#include "ivymike/smart_ptr.h"
#include "ivymike/multiple_alignment.h"

using namespace ivy_mike;
using namespace ivy_mike::tree_parser_ms;

using namespace boost::numeric;

using sequence_model::tag_aa;
using sequence_model::tag_dna;
using sequence_model::model;

template<typename seq_model>
class vu_config {
};


template<>
class vu_config<tag_dna> {
public:
    const static size_t width = 8;
    typedef short scalar;
};

template<>
class vu_config<tag_aa> {
public:
    const static size_t width = 4;
    typedef int scalar;
};

namespace {

const int score_gap_open = 3;
const int score_gap_extend = 1;
const int score_mismatch = 1;
const int score_match_cgap = 3;


//typedef sequence_model::model<sequence_model::tag_dna> seq_model;


//uint8_t normalize_dna( uint8_t c ) {
////    c = std::toupper(c);
////
////    if( c == 'U' ) {
////        c = 'T';
////    }
////
////    return c;
//
//    return seq_model::normalize(c);
//}

char num_to_ascii( int n ) {
    if( n >= 0 && n <= 9 ) {
        return '0' + n;
    } else if( n >= 0xa && n <= 0xf ) {
        return 'a' + n;
    } else {
        throw std::runtime_error( "not a single digit (hex) number" );
    }
}


}


namespace {
    typedef boost::iostreams::tee_device<std::ostream, std::ofstream> log_device;
    typedef boost::iostreams::stream<log_device> log_stream;
    
    log_stream lout;
    
    template<typename stream_, typename device_>
    class bios_open_guard {
        stream_ &m_stream;
    public:
        bios_open_guard( stream_ &stream, device_ &device ) : m_stream(stream) {
            m_stream.open( device );
        }
        ~bios_open_guard() {
            m_stream.close();
        }
    };
    
    typedef bios_open_guard<log_stream, log_device> log_stream_guard;
}

class ostream_test {
    std::ostream &m_os;

public:
    ostream_test( std::ostream &os ) : m_os(os) {}
    void operator()(int i) {
        m_os << i;
    }
};



static bool g_dump_aux = false;



template<class pvec_t,typename seq_tag>
class my_adata_gen : public ivy_mike::tree_parser_ms::adata {
//     static int ct;
    //std::vector<parsimony_state> m_pvec;
    pvec_t m_pvec;

    //typedef seq_model seq_model_t;

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


        m_pvec.init2( seq, model<seq_tag>() );
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



uint8_t to_hex( double v ) {
    int vi = int(fabs(v));

    vi = std::min( vi, 15 );

    if( vi <= 9 ) {
        return '0' + vi;
    } else {
        return 'a' + (vi - 10);
    }

}


template<class pvec_t, typename seq_tag>
void do_newview( pvec_t &root_pvec, lnode *n1, lnode *n2, bool incremental ) {
    typedef my_adata_gen<pvec_t, seq_tag > my_adata;

    std::deque<rooted_bifurcation<lnode> > trav_order;

    //std::cout << "traversal for branch: " << *(n1->m_data) << " " << *(n2->m_data) << "\n";

    rooted_traversal_order( n1, n2, trav_order, incremental );
//     std::cout << "traversal: " << trav_order.size() << "\n";

    for( std::deque< rooted_bifurcation< ivy_mike::tree_parser_ms::lnode > >::iterator it = trav_order.begin(); it != trav_order.end(); ++it ) {
//         std::cout << *it << "\n";

        my_adata *p = it->parent->m_data->get_as<my_adata>();
        my_adata *c1 = it->child1->m_data->get_as<my_adata>();
        my_adata *c2 = it->child2->m_data->get_as<my_adata>();
//         rooted_bifurcation<ivy_mike::tree_parser_ms::lnode>::tip_case tc = it->tc;



//         std::cout << "tip case: " << (*it) << "\n";
        pvec_t::newview(p->get_pvec(), c1->get_pvec(), c2->get_pvec(), it->child1->backLen, it->child2->backLen, it->tc);

    }





    {
        my_adata *c1 = n1->m_data->get_as<my_adata>();
        my_adata *c2 = n2->m_data->get_as<my_adata>();

//         tip_case tc;

        if( c1->isTip && c2->isTip ) {
//                 std::cout << "root: TIP TIP\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_TIP );
        } else if( c1->isTip && !c2->isTip ) {
//                 std::cout << "root: TIP INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c2->get_pvec();
        } else if( !c1->isTip && c2->isTip ) {
//                 std::cout << "root: INNER TIP\n";
            pvec_t::newview(root_pvec, c2->get_pvec(), c1->get_pvec(), n1->backLen, n2->backLen, TIP_INNER );
//             root_pvec = c1->get_pvec();
        } else {
//                 std::cout << "root: INNER INNER\n";
            pvec_t::newview(root_pvec, c1->get_pvec(), c2->get_pvec(), n1->backLen, n2->backLen, INNER_INNER );
        }


    }
//     std::cout << std::hex;
//     for( std::vector< parsimony_state >::const_iterator it = root_pvec.begin(); it != root_pvec.end(); ++it ) {
//         std::cout << *it;
//     }
//
//     std::cout << std::dec << std::endl;

}


//static void seq_to_nongappy_pvec( std::vector<uint8_t> &seq, std::vector<uint8_t> &pvec ) {
//    pvec.resize( 0 );
//
//    for( unsigned int i = 0; i < seq.size(); i++ ) {
//        seq_model::pars_state_t ps = seq_model::s2p(seq[i]);
//
//        if( seq_model::is_single(ps)) {
//            pvec.push_back(ps);
//        }
//
//    }
//
//}

void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq );


template<typename seq_tag>
class queries {
    typedef model<seq_tag> seq_model;

public:

    queries( const char *opt_qs_name ) {


        if( opt_qs_name != 0 ) {
            //
            // read query sequences
            //

            if( opt_qs_name != 0 ) {
                std::ifstream qsf( opt_qs_name );

                if( qsf.bad() ) {
                    throw std::runtime_error( "cannot open qs file");
                }

                // mix them with the qs from the ref alignment
                read_fasta( qsf, m_qs_names, m_qs_seqs);
            }

            if( m_qs_names.empty() ) {
                throw std::runtime_error( "no qs" );
            }

            //
            // setup qs best-score/best-edge lists
            //


            m_qs_pvecs.resize( m_qs_names.size() );
        }
//        m_qs_bestscore.resize(m_qs_names.size());
//        std::fill( m_qs_bestscore.begin(), m_qs_bestscore.end(), 32000);
//        m_qs_bestedge.resize(m_qs_names.size());




    }


    // WARNING: unsafe move semantics on qs
    void add( const std::string &name, std::vector<uint8_t> *qs ) {
        m_qs_names.push_back(name);
        ivy_mike::push_back_swap(m_qs_seqs, *qs );
    }

    void preprocess() {
        //
        // preprocess query sequences
        //
        if( m_qs_seqs.empty() ) {
            throw std::runtime_error( "no query sequences" );
        }

        assert( m_qs_seqs.size() == m_qs_names.size() );
        m_qs_pvecs.resize(m_qs_seqs.size());

        for( size_t i = 0; i < m_qs_seqs.size(); i++ ) {
//            seq_to_nongappy_pvec( m_qs_seqs[i], m_qs_pvecs[i] );
  //          static void seq_to_nongappy_pvec( std::vector<uint8_t> &seq, std::vector<uint8_t> &pvec ) {

            // the following line means: transform sequence to pvec using seq_model::s2p as mapping
            // function and only append the mapped character into pvec, if it corresponds to a single (=non gap)
            // character.

            std::transform( m_qs_seqs[i].begin(), m_qs_seqs[i].end(),
                    back_insert_ifer(m_qs_pvecs[i], seq_model::is_single),
                    seq_model::s2p );

//            for( unsigned int i = 0; i < seq.size(); i++ ) {
//                seq_model::pars_state_t ps = seq_model::s2p(seq[i]);
//
//                if( seq_model::is_single(ps)) {
//                    pvec.push_back(ps);
//                }
//
//            }



        }

//        if( write_testbench ) {
//
//            write_qs_pvecs( "qs.bin" );
//            write_ref_pvecs( "ref.bin" );
//        }

    }


    size_t size() const {
        return m_qs_seqs.size();
    }

    const std::string &name_at( size_t i ) const {
        return m_qs_names.at(i);
    }

    const std::vector<uint8_t> &pvec_at( size_t i ) const {
        assert( m_qs_pvecs.size() == m_qs_seqs.size() );

        return m_qs_pvecs.at(i);
    }

    const std::vector<uint8_t> &seq_at( size_t i ) const {
        return m_qs_seqs.at(i);
    }

    void write_pvecs( const char * name ) {
        std::ofstream os( name );

        os << m_qs_pvecs.size();
        for( std::vector< std::vector< uint8_t > >::iterator it = m_qs_pvecs.begin(); it != m_qs_pvecs.end(); ++it ) {
            os << " " << it->size() << " ";
            os.write( (char *)it->data(), it->size() );

        }
    }


    size_t max_name_length() const {
        size_t len = 0;
        for( std::vector <std::string >::const_iterator it = m_qs_names.begin(); it != m_qs_names.end(); ++it ) {
            len = std::max( len, it->size() );
        }

        return len;
    }


private:
    std::vector <std::string> m_qs_names;
    std::vector <std::vector<uint8_t> > m_qs_seqs;

    std::vector<std::vector <uint8_t> > m_qs_pvecs;

};


template<typename pvec_t, typename seq_tag>
class references {
public:

    typedef model<seq_tag> seq_model;
    typedef my_adata_gen<pvec_t,seq_tag> my_adata;


    references( const char* opt_tree_name, const char *opt_alignment_name, queries<seq_tag> *qs )
      : m_ln_pool(new ln_pool( std::auto_ptr<node_data_factory>(new my_fact<my_adata>) )), spg_(pvec_pgap::pgap_model, &pm_)
    {

            //std::cerr << "papara_nt instantiated as: " << typeid(*this).name() << "\n";
        lout << "papara_nt instantiated as: " << ivy_mike::demangle(typeid(*this).name()) << "\n";
        lout << "scores: " << score_gap_open << " " << score_gap_extend << " " << score_mismatch << " " << score_match_cgap << "\n";



//        std::cerr << ivy_mike::isa<papara_nt<pvec_cgap> >(*this) << " " << ivy_mike::isa<papara_nt<pvec_pgap> >(*this) << "\n";
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
//             std::cout << (*it)->m_data->tipName << "\n";
            name_to_lnode[(*it)->m_data->tipName] = *it;
        }


        {
            //
            // read reference alignment: store the ref-seqs in the tips of the ref-tree
            //
            multiple_alignment ref_ma;
            ref_ma.load_phylip( opt_alignment_name );



            for( unsigned int i = 0; i < ref_ma.names.size(); i++ ) {

                std::map< std::string, sptr::shared_ptr<lnode> >::iterator it = name_to_lnode.find(ref_ma.names[i]);

                // process sequences fomr the ref_ma depending on, if they are contained in the tree.
                // if they are, they are 'swapped' into m_ref_seqs
                // if they are not, into m_qs_seqs. (gaps in the QS are removed later)

                if( it != name_to_lnode.end() ) {
                    sptr::shared_ptr< lnode > ln = it->second;
                    //      adata *ad = ln->m_data.get();

                    assert( ivy_mike::isa<my_adata>(ln->m_data.get()) ); //typeid(*ln->m_data.get()) == typeid(my_adata ) );
                    my_adata *adata = static_cast<my_adata *> (ln->m_data.get());

                    m_ref_names.push_back(std::string() );
                    m_ref_seqs.push_back(std::vector<uint8_t>() );

                    m_ref_names.back().swap( ref_ma.names[i] );
                    m_ref_seqs.back().swap( ref_ma.data[i] );

                    // WARNING: make sure not to keep references to elements of m_ref_seqs at this point!
                    adata->init_pvec( m_ref_seqs.back() );
                } else {
                    qs->add(ref_ma.names[i], &ref_ma.data[i]);
                }
            }
        }
        pm_.reset( m_ref_seqs );
        std::cout << "p: " << pm_.setup_pmatrix(0.1) << "\n";

        //
        // collect list of edges
        //

        visit_edges( n, m_ec );

        lout << "edges: " << m_ec.m_edges.size() << "\n";

    }

    void build_ref_vecs() {
        // pre-create the ancestral state vectors. This step is necessary for the threaded version, because otherwise, each
        // thread would need an independent copy of the tree to do concurrent newviews. Anyway, having a copy of the tree
        // in each thread will most likely use more memory than storing the pre-calculated vectors.

        // TODO: maybe try lazy create/cache of the asv's in the threads

        ivy_mike::timer t1;



        assert( m_ref_aux.empty() && m_ref_pvecs.empty() );

        m_ref_pvecs.reserve( m_ec.m_edges.size() );
        m_ref_aux.reserve( m_ec.m_edges.size() );


        for( size_t i = 0; i < m_ec.m_edges.size(); i++ ) {
            pvec_t root_pvec;

//             std::cout << "newview for branch " << i << ": " << *(m_ec.m_edges[i].first->m_data) << " " << *(m_ec.m_edges[i].second->m_data) << "\n";

            if( i == 340 ) {
                g_dump_aux = true;
            }

            do_newview<pvec_t,seq_tag>( root_pvec, m_ec.m_edges[i].first, m_ec.m_edges[i].second, true );

            g_dump_aux = false;
            // TODO: try something fancy with rvalue refs...

            m_ref_pvecs.push_back( std::vector<int>() );
            m_ref_aux.push_back( std::vector<unsigned int>() );

            root_pvec.to_int_vec(m_ref_pvecs.back());
            root_pvec.to_aux_vec(m_ref_aux.back());


            m_ref_gapp.push_back( std::vector<double>() );

            if( ivy_mike::same_type<pvec_t,pvec_pgap>::result ) {
                // WTF: this is why mixing static and dynamic polymorphism is a BAD idea!
                pvec_pgap *rvp = reinterpret_cast<pvec_pgap *>(&root_pvec);
                rvp->to_gap_post_vec(m_ref_gapp.back());

//              std::transform( m_ref_gapp.back().begin(), m_ref_gapp.back().end(), std::ostream_iterator<int>(std::cout), ivy_mike::scaler_clamp<double>(10,0,9) );
//
//              std::cout << "\n";
            }


        }

        std::cout << "pvecs created: " << t1.elapsed() << "\n";

    }


    const std::vector<int> &pvec_at( size_t i ) const {
        return m_ref_pvecs.at(i);
    }

    const std::vector<unsigned int> &aux_at( size_t i ) const {
        return m_ref_aux.at(i);
    }

    size_t num_pvecs() const {
        return m_ref_pvecs.size();
    }

    size_t pvec_size() const {
        assert( !m_ref_pvecs.empty());
        return m_ref_pvecs.front().size();
    }

    void write_pvecs( const char * name ) {
        std::ofstream os( name );

        os << m_ref_pvecs.size();
        for( size_t i = 0; i < m_ref_pvecs.size(); ++i ) {
            os << " " << m_ref_pvecs[i].size() << " ";
            os.write( (char *)m_ref_pvecs[i].data(), m_ref_pvecs[i].size() * sizeof(int));
            os.write( (char *)m_ref_aux[i].data(), m_ref_aux[i].size() * sizeof(unsigned int));
        }
    }


    size_t max_name_length() const {
        size_t len = 0;
        for( std::vector <std::string >::const_iterator it = m_ref_names.begin(); it != m_ref_names.end(); ++it ) {
            len = std::max( len, it->size() );
        }

        return len;
    }


    void write_seqs( std::ostream &os, size_t pad ) {
        for( size_t i = 0; i < m_ref_seqs.size(); i++ ) {
            os << std::setw(pad) << std::left << m_ref_names[i];
            std::transform( m_ref_seqs[i].begin(), m_ref_seqs[i].end(), std::ostream_iterator<char>(os), seq_model::normalize );
            os << "\n";
        }
    }
private:
    std::vector <std::string > m_ref_names;
    std::vector <std::vector<uint8_t> > m_ref_seqs;
    std::auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    edge_collector<lnode> m_ec;

    std::vector<std::vector <int> > m_ref_pvecs;
    std::vector<std::vector <unsigned int> > m_ref_aux;
    std::vector<std::vector <double> > m_ref_gapp;

    probgap_model pm_;
    stupid_ptr_guard<probgap_model> spg_;

};


template<typename seq_tag>
class block_queue {
    const static size_t VW = vu_config<seq_tag>::width;

public:
    struct block_t {
        block_t() {
            memset( this, 0, sizeof( block_t )); // FIXME: hmm, this is still legal?
        }

        // WARNING: these are pointers into m_ref_pvecs and m_ref_aux
        // make sure they stay valid!
        const int *seqptrs[VW];
        const unsigned int *auxptrs[VW];
        const double *gapp_ptrs[VW];
        size_t ref_len;
        size_t edges[VW];
        int num_valid;
    };


//    bool empty() {
//        ivy_mike::lock_guard<ivy_mike::mutex> lock(m_qmtx);
//
//        return m_blockqueue.empty();
//
//    }

    bool get_block( block_t *block ) {
        ivy_mike::lock_guard<ivy_mike::mutex> lock( m_qmtx );

        if( m_blockqueue.empty() ) {
            return false;
        }

        *block = m_blockqueue.front();
        m_blockqueue.pop_front();

        return true;

    }


    // WARNING: this method is not synchronized, and shall only be called before the worker threads are running
    void push_back( const block_t &b ) {
        m_blockqueue.push_back(b);
    }

    ivy_mike::mutex *hack_mutex() {
        return &m_qmtx;
    }
private:
    ivy_mike::mutex m_qmtx; // mutex for the block queue and the qs best score/edge arrays
    std::deque<block_t> m_blockqueue;
    std::vector <int> m_qs_bestscore;
    std::vector <int> m_qs_bestedge;
};


class scoring_results {

public:
    scoring_results( size_t num_qs )
    : best_score_(num_qs, std::numeric_limits<int>::max() ),
      best_ref_(num_qs, size_t(-1))
    {}



    bool offer( size_t qs, size_t ref, int score ) {
        ivy_mike::lock_guard<ivy_mike::mutex> lock(mtx_);

        if( best_score_.at(qs) > score || (best_score_.at(qs) == score && ref < best_ref_.at(qs))) {
            best_score_[qs] = score;
            best_ref_.at(qs) = ref;
            return true;
        }

        return false;

    }


    template<typename idx_iter, typename score_iter>
    void offer( size_t qs, idx_iter ref_start, idx_iter ref_end, score_iter score_start ) {
        ivy_mike::lock_guard<ivy_mike::mutex> lock(mtx_);


        while( ref_start != ref_end ) {
            if( best_score_.at(qs) > *score_start || (best_score_.at(qs) == *score_start && *ref_start < best_ref_.at(qs))) {
                best_score_[qs] = *score_start;
                best_ref_.at(qs) = *ref_start;
            }


            ++ref_start;
            ++score_start;
        }

    }


    int bestscore_at(size_t i ) const {
        return best_score_.at(i);
    }

    size_t bestedge_at(size_t i ) const {
        return best_ref_.at(i);
    }

private:
    std::vector<int> best_score_;
    std::vector<size_t> best_ref_;

    ivy_mike::mutex mtx_;

};


template<typename seq_tag>
class worker {



    const static size_t VW = vu_config<seq_tag>::width;
    typedef typename vu_config<seq_tag>::scalar vu_scalar_t;
    typedef typename block_queue<seq_tag>::block_t block_t;




    block_queue<seq_tag> &block_queue_;
    scoring_results &results_;

    const queries<seq_tag> &qs_;

    size_t m_rank;
public:
    worker( block_queue<seq_tag> *bq, scoring_results *res, const queries<seq_tag> &qs, size_t rank ) : block_queue_(*bq), results_(*res), qs_(qs), m_rank(rank) {}
    void operator()() {

        pars_align_seq::arrays seq_arrays(true);
        pars_align_gapp_seq::arrays seq_arrays_gapp(true);

        ivy_mike::timer tstatus;
        double last_tstatus = 0.0;

        uint64_t ncup = 0;
        while( true ) {
            block_t block;

            if( !block_queue_.get_block(&block)) {
                break;
            }
#if 1
       //     assert( VW == 8 );
            const align_pvec_score<vu_scalar_t,VW> aligner( block.seqptrs, block.auxptrs, block.ref_len, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
            for( unsigned int i = 0; i < qs_.size(); i++ ) {


                aligner.align(qs_.pvec_at(i));
                const vu_scalar_t *score_vec = aligner.get_scores();

                ncup += block.num_valid * block.ref_len * qs_.pvec_at(i).size();

                results_.offer( i, block.edges, block.edges + block.num_valid, score_vec );

            }
        }
#else

//        assert( block.gapp_ptrs[0] != 0 );
//        assert( VW == 4 );
//        const align_pvec_gapp_score<4> aligner( block.seqptrs, block.gapp_ptrs, block.ref_len, score_mismatch, score_match_cgap, score_gap_open, score_gap_extend );
//        for( unsigned int i = 0; i < m_pnt.m_qs_names.size(); i++ ) {
//
//            size_t stride = 1;
//            size_t aux_stride = 1;
//
//            aligner.align(m_pnt.m_qs_pvecs[i]);
//            const float *score_vec = aligner.get_scores();
//
//            ncup += block.num_valid * block.ref_len * m_pnt.m_qs_pvecs[i].size();
//            {
//                ivy_mike::lock_guard<ivy_mike::mutex> lock( m_pnt.m_qmtx );
//
//                for( int k = 0; k < block.num_valid; k++ ) {
//
//
//
//                    if( score_vec[k] < m_pnt.m_qs_bestscore[i] || (score_vec[k] == m_pnt.m_qs_bestscore[i] && block.edges[k] < m_pnt.m_qs_bestedge[i] )) {
//                        const bool validate = false;
//                        if( validate ) {
//                            const int *seqptr = block.seqptrs[k];
//                            const double *gapp_ptr = block.gapp_ptrs[k];
//
////                                std::vector<double> gapp_tmp(gapp_ptr, gapp_ptr + block.ref_len);
//
//
//                            pars_align_gapp_seq pas( seqptr, m_pnt.m_qs_pvecs[i].data(), block.ref_len, m_pnt.m_qs_pvecs[i].size(), stride, gapp_ptr, aux_stride, seq_arrays_gapp, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
//                            int res = pas.alignFreeshift(INT_MAX);
//
//                            if( res != score_vec[k] ) {
//
//
//                                std::cout << "meeeeeeep! score: " << score_vec[k] << " " << res << "\n";
//                            }
//                        }
//
//                        m_pnt.m_qs_bestscore[i] = score_vec[k];
//                        m_pnt.m_qs_bestedge[i] = block.edges[k];
//                    }
//                }
//            }
//        }
//    }

#endif
        {
            ivy_mike::lock_guard<ivy_mike::mutex> lock( *block_queue_.hack_mutex() );
            std::cout << "thread " << m_rank << ": " << ncup / (tstatus.elapsed() * 1e9) << " gncup/s\n";
        }
    }
};


template<typename pvec_t,typename seq_tag>
void build_block_queue( const references<pvec_t,seq_tag> &refs, block_queue<seq_tag> *bq ) {
    // creates the list of ref-block to be consumed by the worker threads.  A ref-block onsists of N ancestral state sequences, where N='width of the vector unit'.
    // The vectorized alignment implementation will align a QS against a whole ref-block at a time, rather than a single ancestral state sequence as in the
    // sequencial algorithm.

    const static size_t VW = vu_config<seq_tag>::width;

    typedef typename block_queue<seq_tag>::block_t block_t;


    size_t n_groups = (refs.num_pvecs() / VW);
    if( (refs.num_pvecs() % VW) != 0 ) {
        n_groups++;
    }


//         std::vector<int> seqlist[VW];
//         const int *seqptrs[VW];
//         std::vector<unsigned int> auxlist[VW];
//         const unsigned int *auxptrs[VW];



    for ( size_t j = 0; j < n_groups; j++ ) {
        int num_valid = 0;



        block_t block;

        for( unsigned int i = 0; i < VW; i++ ) {

            size_t edge = j * VW + i;
            if( edge < refs.num_pvecs()) {
                block.edges[i] = edge;
                block.num_valid++;

                block.seqptrs[i] = refs.pvec_at(edge).data();
                block.auxptrs[i] = refs.aux_at(edge).data();

//                if( !m_ref_gapp[edge].empty() ) {
//                    block.gapp_ptrs[i] = m_ref_gapp[edge].data();
//                } else {
                block.gapp_ptrs[i] = 0;
//                }

                block.ref_len = refs.pvec_size();
                //                     do_newview( root_pvec, m_ec.m_edges[edge].first, m_ec.m_edges[edge].second, true );
//                     root_pvec.to_int_vec(seqlist[i]);
//                     root_pvec.to_aux_vec(auxlist[i]);
//
//                     seqptrs[i] = seqlist[i].data();
//                     auxptrs[i] = auxlist[i].data();

                num_valid++;
            } else {
                if( i < 1 ) {
                    std::cout << "edge: " << edge << " " << refs.num_pvecs() << std::endl;

                    throw std::runtime_error( "bad integer mathematics" );
                }
                block.edges[i] = block.edges[i-1];

                block.seqptrs[i] = block.seqptrs[i-1];
                block.auxptrs[i] = block.auxptrs[i-1];
                block.gapp_ptrs[i] = block.gapp_ptrs[i-1];
            }

        }
        bq->push_back(block);
    }
}


template<typename seq_tag>
void seq_to_position_map(const std::vector< uint8_t >& seq, std::vector< int > &map) {
    typedef model<seq_tag> seq_model;

    for( size_t i = 0; i < seq.size(); ++i ) {
        if( seq_model::is_single(seq_model::s2p(seq[i]))) {
            map.push_back(int(i));
        }
    }
}

void gapstream_to_position_map( const std::vector< uint8_t >& gaps, std::vector< int > &map) {
    align_utils::trace_to_position_map( gaps, &map );

}




template<typename pvec_t, typename seq_tag>
void calc_scores( size_t n_threads, const references<pvec_t, seq_tag> &refs, const queries<seq_tag> &qs, scoring_results *res ) {

    //
    // build the alignment blocks
    //


    block_queue<seq_tag> bq;
    build_block_queue(refs, &bq);

    //
    // work
    //
    ivy_mike::timer t1;
    ivy_mike::thread_group tg;
    lout << "start scoring, using " << n_threads <<  " threads" << std::endl;

    typedef worker<seq_tag> worker_t;

    for( size_t i = 1; i < n_threads; ++i ) {
        tg.create_thread(worker_t(&bq, res, qs, i));
    }

    worker_t w0(&bq, res, qs, 0 );
    w0();

    tg.join_all();

    lout << "scoring finished: " << t1.elapsed() << std::endl;

}
template<typename seq_tag>
void print_best_scores( std::ostream &os, const queries<seq_tag> &qs, const scoring_results &res ) {
    boost::io::ios_all_saver ioss(os);
    os << std::setfill ('0');
    for( unsigned int i = 0; i < qs.size(); i++ ) {
        os << qs.name_at(i) << " "  << std::setw (4) << res.bestedge_at(i) << " " << std::setw(5) << res.bestscore_at(i) << "\n";
    }
}


void gapstream_to_alignment( const std::vector<uint8_t> &gaps, const std::vector<uint8_t> &raw, std::vector<uint8_t> *out, uint8_t gap_char ) {

    std::vector<uint8_t>::const_reverse_iterator rit = raw.rbegin();

    for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {

        if ( *git == 1) {
            out->push_back(gap_char);
        } else if ( *git == 0 ) {
            assert( rit < raw.rend() );
            out->push_back(*rit);
            ++rit;
        } else {
            ++rit; // just consume one QS character
        }
    }

    std::reverse(out->begin(), out->end());
}



template<typename pvec_t, typename seq_tag>
void align_best_scores( std::ostream &os, std::ostream &os_quality, const queries<seq_tag> &qs, const references<pvec_t,seq_tag> &refs, const scoring_results &res ) {
    // create the actual alignments for the best scoring insertion position (=do the traceback)

    lout << "generating best scoring alignments\n";
    ivy_mike::timer t1;
    pars_align_seq::arrays seq_arrays(true);
    pars_align_gapp_seq::arrays seq_arrays_gapp(true);

    double mean_quality = 0.0;
    double n_quality = 0.0;

    std::vector<uint8_t> out_qs;

    for( size_t i = 0; i < qs.size(); i++ ) {
        int best_edge = res.bestedge_at(i);

        assert( best_edge >= 0 && size_t(best_edge) < refs.num_pvecs() );

        int score = -1;
        std::vector<uint8_t> tbv;

//        if( true ) {

        const std::vector<uint8_t> &qp = qs.pvec_at(i);
        const int *seqptr = refs.pvec_at(best_edge).data();
        const unsigned int *auxptr = refs.aux_at(best_edge).data();

        const size_t ref_len = refs.pvec_size();

        const size_t stride = 1;
        const size_t aux_stride = 1;
        pars_align_seq pas( seqptr, qp.data(), ref_len, qp.size(), stride, auxptr, aux_stride, seq_arrays, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
        score = pas.alignFreeshift(INT_MAX);
        pas.tracebackCompressed(tbv);
//        } else {
//            const int *seqptr = m_ref_pvecs[best_edge].data();
//
//            assert( !m_ref_gapp[best_edge].empty() );
//
//            const double *gappptr = m_ref_gapp[best_edge].data();
//
//            const size_t ref_len = m_ref_pvecs[best_edge].size();
//
//            const size_t stride = 1;
//            const size_t aux_stride = 1;
//            pars_align_gapp_seq pas( seqptr, m_qs_pvecs[i].data(), ref_len, m_qs_pvecs[i].size(), stride, gappptr, aux_stride, seq_arrays_gapp, 0, score_gap_open, score_gap_extend, score_mismatch, score_match_cgap );
//            res = pas.alignFreeshift(INT_MAX);
//            pas.tracebackCompressed(tbv);
//        }
//            std::copy( tbv.begin(), tbv.end(), std::ostream_iterator<int>(std::cout, " " ));
//            std::cout << "\n";





        out_qs.clear();

        gapstream_to_alignment(tbv, qp, &out_qs, 0xf);

        //
        // in place transform out_qs from pstate to dna form. WATCH OUT!
        //

        std::transform( out_qs.begin(), out_qs.end(), out_qs.begin(), dna_parsimony_mapping::p2d );

        os << qs.name_at(i) << "\t";
        std::copy( out_qs.begin(), out_qs.end(), std::ostream_iterator<char>(os));
        os << "\n";



        if( score != res.bestscore_at(i) ) {
            std::cout << "meeeeeeep! score: " << res.bestscore_at(i) << " " << score << "\n";
        }

        if( os_quality.good() && out_qs.size() == qs.seq_at(i).size() ) {


            std::vector<int> map_ref;
            std::vector<int> map_aligned;
            seq_to_position_map<seq_tag>( qs.seq_at(i), map_ref );
            align_utils::trace_to_position_map( tbv, &map_aligned );


            if( map_ref.size() != map_aligned.size() ) {
                throw std::runtime_error( "alignment quirk: map_ref.size() != map_aligned.size()" );
            }

            size_t num_equal = ivy_mike::count_equal( map_ref.begin(), map_ref.end(), map_aligned.begin() );

            //std::cout << "size: " << map_ref.size() << " " << map_aligned.size() << " " << m_qs_seqs[i].size() << "\n";
            //std::cout << num_equal << " equal of " << map_ref.size() << "\n";

            double score = num_equal / double(map_ref.size());
            //double score = alignment_quality( out_qs, m_qs_seqs[i], debug );

            os_quality << qs.name_at(i) << " " << score << "\n";

            mean_quality += score;
            n_quality += 1;
        }


    }
    lout << "alignment finished: " << t1.elapsed() << "\n";
    lout << "mean quality: " << mean_quality / n_quality << "\n";

}


double alignment_quality_very_strict ( const std::vector< uint8_t > &s1, const std::vector< uint8_t >& s2, bool debug = false ) {
    size_t nident = 0;
    size_t ngap1 = 0;
    size_t ngap2 = 0;


    for( std::vector< uint8_t >::const_iterator it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2 ) {

        if( dna_parsimony_mapping::is_gap( *it1 ) ) {
            ngap1++;
        }

        if( dna_parsimony_mapping::is_gap( *it2 ) ) {
            ngap2++;
        }
        if( debug ) {
            std::cerr << ngap1 << " " << ngap2 << " " << *it1 << " " << *it2 << "\n";
        }

        if( ngap1 == ngap2 ) {
            nident++;
        }
    }

    return double(nident) / s1.size();

}
double alignment_quality ( const std::vector< uint8_t > &s1, const std::vector< uint8_t >& s2, bool debug = false ) {
    size_t nident = 0;

//         size_t nident_nongap = 0;
//         size_t n_nongap = 0;

    for( std::vector< uint8_t >::const_iterator it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2 ) {
        if( dna_parsimony_mapping::d2p(*it1) == dna_parsimony_mapping::d2p(*it2) ) {
            nident++;
        }
    }

    return double(nident) / s1.size();

}



std::string filename( const std::string &run_name, const char *type ) {
    std::stringstream ss;
    
    ss << "papara_" << type << "." << run_name;
    
    return ss.str();
}

bool file_exists(const char *filename)
{
  std::ifstream is(filename);
  return is;
}

template<typename pvec_t, typename seq_tag>
void run_papara( const std::string &qs_name, const std::string &alignment_name, const std::string &tree_name, size_t num_threads, const std::string &run_name ) {

    queries<seq_tag> qs(qs_name.c_str());
    references<pvec_t,seq_tag> refs( tree_name.c_str(), alignment_name.c_str(), &qs );

    qs.preprocess();
    refs.build_ref_vecs();

    scoring_results res( qs.size() );


    const size_t VW = 8;

    calc_scores<pvec_t, seq_tag>(num_threads, refs, qs, &res );

    std::string score_file(filename(run_name, "alignments"));
    std::string quality_file(filename(run_name, "quality"));


    size_t pad = 1 + std::max(qs.max_name_length(), refs.max_name_length());

    std::ofstream os( score_file.c_str() );
    assert( os.good() );

    std::ofstream os_qual( quality_file.c_str() );
    assert( os_qual.good() );

    refs.write_seqs(os, pad);
    align_best_scores( os, os_qual, qs, refs, res );

}

int main( int argc, char *argv[] ) {

//     aligned_buffer<int> xxx(1024);
    
    

    
    namespace igo = ivy_mike::getopt;

    ivy_mike::getopt::parser igp;

    std::string opt_tree_name;
    std::string opt_alignment_name;
    std::string opt_qs_name;
    bool opt_use_cgap;
    int opt_num_threads;
    std::string opt_run_name;
    bool opt_write_testbench;

    
    igp.add_opt( 't', igo::value<std::string>(opt_tree_name) );
    igp.add_opt( 's', igo::value<std::string>(opt_alignment_name) );
    igp.add_opt( 'q', igo::value<std::string>(opt_qs_name) );
    igp.add_opt( 'c', igo::value<bool>(opt_use_cgap, true).set_default(false) );
    igp.add_opt( 'j', igo::value<int>(opt_num_threads).set_default(1) );
    igp.add_opt( 'n', igo::value<std::string>(opt_run_name).set_default("default") );
    igp.add_opt( 'b', igo::value<bool>(opt_write_testbench, true).set_default(false) );
    
    igp.parse(argc,argv);

    if( igp.opt_count('t') != 1 || igp.opt_count('s') != 1  ) {
        std::cerr << "missing options -t and/or -s (-q is optional)\n";
        return 0;
    }
    ivy_mike::timer t;

    const char *qs_name = 0;
    if( !opt_qs_name.empty() ) {
        qs_name = opt_qs_name.c_str();
    }
    
    std::string log_filename = filename( opt_run_name, "log" );

    if( opt_run_name != "default" && file_exists(log_filename.c_str()) ) {
        std::cout << "log file already exists for run '" << opt_run_name << "'\n";
        return 0;
    }
    
    std::ofstream logs( log_filename.c_str());
    if( !logs ) {
        std::cout << "could not open logfile for writing: " << log_filename << std::endl;
        return 0;
    }
    
    log_device ldev( std::cout, logs );
    log_stream_guard lout_guard( lout, ldev );
    
    




    if( opt_use_cgap ) {

        run_papara<pvec_pgap, tag_aa>( qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name );
    } else {
        run_papara<pvec_pgap, tag_dna>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name );
    }

    
    
    std::cout << t.elapsed() << std::endl;
    lout << "SUCCESS " << t.elapsed() << std::endl;
    return 0;
//     getchar();
}

