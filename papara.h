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

#ifndef __papara_h
#define __papara_h

#include <ivymike/disable_shit.h>
#define BOOST_UBLAS_NDEBUG 1

// #include <stdexcept>
#include <iostream>
//
#include <vector>
#include <deque>
#include <map>
#include <functional>
#include <cstring>
#include <algorithm>
#include <numeric>

#include <boost/io/ios_state.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>




#include "sequence_model.h"
// #include "parsimony.h"
#include "pvec.h"
// #include "align_utils.h"




// #include "align_utils.h"

#include "ivymike/tree_parser.h"

#include "ivymike/thread.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
#include "ivymike/smart_ptr.h"
#include "ivymike/multiple_alignment.h"


namespace papara {
using namespace ivy_mike;
// using namespace ivy_mike::tree_parser_ms;

namespace im_tree_parser = ivy_mike::tree_parser_ms;

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
    const static scalar full_mask = scalar(-1);
};

template<>
class vu_config<tag_aa> {
public:
    const static size_t width = 4;
    typedef int scalar;
    const static scalar full_mask = scalar(-1);
};

struct papara_score_parameters {
    
    static papara_score_parameters default_scores() {
        return papara_score_parameters(-3,-1,2,-3);
    }

    static papara_score_parameters parse_scores( const char *opts) {
        
        int o, e, m, mc;
#ifndef _MSC_VER
        int n = sscanf( opts, "%d:%d:%d:%d", &o, &e, &m, &mc );
#else
        int n = sscanf_s( opts, "%d:%d:%d:%d", &o, &e, &m, &mc );
#endif

        if( n != 4 ) {
            std::cerr <<  "cannot parse user options for papara: '" << opts << "'\nIt should match the following format: <open>:<extend>:<match>:<match cg>\n";
            throw std::runtime_error( "bailing out" );
        }
        
        return papara_score_parameters(o, e, m, mc);

    }
    
    papara_score_parameters( int open_, int ext_, int match_, int match_cg_ ) 
     : gap_open( open_ ),
       gap_extend( ext_ ),
       match( match_ ),
       match_cgap( match_cg_ )
    {}
    
    void print( std::ostream &os ) const {
        os << "scoring parameters: " << gap_open << " " << gap_extend << " " << match << " " << match_cgap << "\n";
    }
    
    int gap_open;
    int gap_extend;
    int match;
    int match_cgap;
    
    bool operator==(const papara_score_parameters &other ) {
        return gap_open == other.gap_open && gap_extend == other.gap_extend && match == other.match && match_cgap == other.match_cgap;
    }
    
    bool operator!=( const papara_score_parameters &other ) {
        return !operator==(other);
    }
    
};

namespace {

// const int score_gap_open = -3;
// const int score_gap_extend = -1;
// const int score_match = 2;
// const int score_match_cgap = -3;


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

extern log_stream lout;

class ostream_test {
    std::ostream &m_os;

public:
    ostream_test( std::ostream &os ) : m_os(os) {}
    void operator()(int i) {
        m_os << i;
    }
};



extern bool g_dump_aux;



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


    static void normalize_name( std::string & str ) ;

public:

    typedef typename seq_model::pars_state_t pars_state_t;

    queries( const std::string &opt_qs_name );



    // WARNING: unsafe move semantics on qs
    void add( const std::string &name, std::vector<uint8_t> *qs ) ;

    void preprocess() ;


    size_t size() const {
        return m_qs_seqs.size();
    }

    const std::string &name_at( size_t i ) const {
        return m_qs_names.at(i);
    }

    const std::vector<pars_state_t> &pvec_at( size_t i ) const {
        assert( m_qs_pvecs.size() == m_qs_seqs.size() );

        return m_qs_pvecs.at(i);
    }

    const std::vector<uint8_t> &seq_at( size_t i ) const {
        return m_qs_seqs.at(i);
    }

    const std::vector<uint8_t> &cseq_at( size_t i ) const {
        return m_qs_cseqs.at(i);
    }

    void write_pvecs( const char * name ) ;


    size_t max_name_length() const ;



    size_t calc_cups_per_ref( size_t ref_len ) const ;


private:
    std::vector <std::string> m_qs_names;
    std::vector <std::vector<uint8_t> > m_qs_seqs;

    std::vector <std::vector<uint8_t> > m_qs_cseqs;

    std::vector<std::vector <pars_state_t> > m_qs_pvecs;

};


template<typename pvec_t, typename seq_tag>
class references {
public:

    typedef model<seq_tag> seq_model;
    typedef my_adata_gen<pvec_t,seq_tag> my_adata;



    references( const char* opt_tree_name, const char *opt_alignment_name, queries<seq_tag> *qs )
      ;

    void remove_full_gaps() {
        
    }
    
    void build_ref_vecs() ;


    const std::string &name_at( size_t i ) const {
        return m_ref_names.at(i);
    }

    const std::vector<uint8_t> & seq_at( size_t i ) const {
        return m_ref_seqs.at(i);
    }

    size_t num_seqs() const {
        return m_ref_seqs.size();
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

    void write_pvecs( const char * name ) ;


    size_t max_name_length() const ;


//    void write_seqs( std::ostream &os, size_t pad ) {
//        for( size_t i = 0; i < m_ref_seqs.size(); i++ ) {
//            os << std::setw(pad) << std::left << m_ref_names[i];
//            std::transform( m_ref_seqs[i].begin(), m_ref_seqs[i].end(), std::ostream_iterator<char>(os), seq_model::normalize );
//            os << "\n";
//        }
//    }
private:
    std::vector <std::string > m_ref_names;
    std::vector <std::vector<uint8_t> > m_ref_seqs;
    std::auto_ptr<ivy_mike::tree_parser_ms::ln_pool> m_ln_pool;
    edge_collector<im_tree_parser::lnode> m_ec;

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
    class candidate {

    public:
        candidate( int score, size_t ref ) : score_(score), ref_(ref) {}

        bool operator<(const candidate &other ) const ;

        int score() const {
            return score_;
        }
        size_t ref() const {
            return ref_;
        }

    private:
        int score_;
        size_t ref_;
    };

    class candidates : private std::vector<candidate> {
    public:

        candidates( size_t max_num )
         : max_num_(max_num)
        {
            reserve(max_num_);
        }

        void offer( int score, size_t ref ) ;

        using std::vector<candidate>::at;
        using std::vector<candidate>::operator[];
        using std::vector<candidate>::size;

    private:
        const size_t max_num_;

        std::vector<candidate> cands_;
    };

public:
    scoring_results( size_t num_qs, const candidates &cands_template )
    : best_score_(num_qs, std::numeric_limits<int>::min() ),
      best_ref_(num_qs, size_t(-1)),
      candss_(num_qs, cands_template )
    {}



    bool offer( size_t qs, size_t ref, int score ) ;


    template<typename idx_iter, typename score_iter>
    void offer( size_t qs, idx_iter ref_start, idx_iter ref_end, score_iter score_start ) {
        ivy_mike::lock_guard<ivy_mike::mutex> lock(mtx_);


        while( ref_start != ref_end ) {
            candss_.at( qs ).offer( *score_start, *ref_start );


            if( best_score_.at(qs) < *score_start || (best_score_.at(qs) == *score_start && *ref_start < best_ref_.at(qs))) {
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

    const candidates &candidates_at( size_t i ) const {
        return candss_.at( i );
    }

private:
    std::vector<int> best_score_;
    std::vector<size_t> best_ref_;

    std::vector<candidates> candss_;

    ivy_mike::mutex mtx_;

};




class ref_gap_collector {
public:

    ref_gap_collector( size_t ref_len ) : ref_gaps_(ref_len + 1) {}

    void add_trace( const std::vector<uint8_t> &gaps ) {

        size_t ptr  = ref_gaps_.size() - 1;

        std::vector<size_t> ref_gaps( ref_gaps_.size() );

        // count how many gaps are inserted before each ref character (the last entry refers to the position after the last ref character)
        for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {


            if( *git == 0 || *git == 1 ) {
                // consume one ref character without inserting gap
                --ptr;
            } else {

                assert( ptr >= 0 );

                // count all gaps inserted at current ref position
                ++ref_gaps[ptr];
            }
        }

        // update the _global_ maximum 'gaps-per-ref-position' map
        std::transform( ref_gaps_.begin(), ref_gaps_.end(), ref_gaps.begin(), ref_gaps_.begin(), std::max<size_t> );
    }

    // TODO: shouldn't it be possible to infer the state_type from oiter?
    template<typename iiter, typename oiter, typename state_type>
    void transform( iiter istart, iiter iend, oiter ostart, state_type gap ) const {

        size_t s = std::distance(istart, iend);
        assert( s == ref_gaps_.size() - 1 );


        size_t i = 0;
        while( istart != iend ) {


            for( size_t j = 0; j < ref_gaps_[i]; ++j ) {
                *(ostart++) = gap;
            }
            *(ostart++) = *(istart++);
            ++i;

        }
        for( size_t j = 0; j < ref_gaps_.back(); ++j ) {
            *(ostart++) = gap;
        }

    }


    size_t gaps_before( size_t i ) const {
        return ref_gaps_.at(i);
    }

    size_t ref_len() const {
        return ref_gaps_.size() - 1;
    }

    size_t transformed_ref_len() const {
        return ref_len() + std::accumulate( ref_gaps_.begin(), ref_gaps_.end(), 0 );
    }

private:

    std::vector<size_t> ref_gaps_;

};


class output_alignment {
public:
    enum seq_type {
        type_ref,
        type_qs
    };
    
    typedef std::vector<char> out_seq;
    
    virtual ~output_alignment() {}
    
    virtual void push_back( const std::string &name, const out_seq &seq, seq_type t ) = 0;
    virtual void set_max_name_length( size_t len ) = 0;
    virtual void set_size( size_t num_rows, size_t num_cols ) = 0;
    
};

class output_alignment_phylip : public output_alignment {
public:
    output_alignment_phylip( const char *filename ) : num_rows_(0), num_cols_(0), max_name_len_(0), header_flushed_(false) {
        os_.open( filename );
        assert( os_.good() );
    }
    
    void set_size( size_t num_rows, size_t num_cols ) {
        num_rows_ = num_rows;
        num_cols_ = num_cols;
    }
    
    void write_seq_phylip( const std::string &name, const out_seq &seq ) ;
    
    void push_back( const std::string &name, const out_seq &seq, seq_type t ) ;
    
    void set_max_name_length( size_t len ) {
        max_name_len_ = len;
    }
private:
    std::ofstream os_;
    size_t num_rows_;
    size_t num_cols_;
    
    size_t max_name_len_; // that's a bad name. already includes the space.
    
    bool header_flushed_;
};


class output_alignment_fasta : public output_alignment {
    
    
public:
    
    
    
    output_alignment_fasta( const char *filename )  {
        os_.open( filename );
        assert( os_.good() );
    }
    
    void set_size( size_t num_rows, size_t num_cols ) {
        
    }
    
    
    
    void push_back( const std::string &name, const out_seq &seq, seq_type t ) ;
    
    void set_max_name_length( size_t len ) {
        
    }
private:
    std::ofstream os_;
    
};


template<typename pvec_t, typename seq_tag>
class driver {
public:
    typedef queries<seq_tag> my_queries;
    typedef references<pvec_t,seq_tag> my_references;
    typedef block_queue<seq_tag> my_block_queue;
    
    static void calc_scores( size_t n_threads, const my_references &refs, const my_queries &qs, scoring_results *res, const papara_score_parameters &sp );
    
    static void do_newview( pvec_t &root_pvec, im_tree_parser::lnode *n1, im_tree_parser::lnode *n2, bool incremental ) ;
    
    static void build_block_queue( const my_references &refs, my_block_queue *bq ) ;
    
    static void seq_to_position_map(const std::vector< uint8_t >& seq, std::vector< int > &map) ;
    
    static void print_best_scores( std::ostream &os, const my_queries &qs, const scoring_results &res ) ;
    
    static std::vector<std::vector<uint8_t> > generate_traces( std::ostream &os_quality, std::ostream &os_cands, const my_queries &qs, const my_references &refs, const scoring_results &res, const papara_score_parameters &sp ) ;
    
    static void align_best_scores2( std::ostream &os, std::ostream &os_quality, std::ostream &os_cands, const my_queries &qs, const my_references &refs, const scoring_results &res, size_t pad, const bool ref_gaps, const papara_score_parameters &sp ) ;
    
    static void align_best_scores( std::ostream &os, std::ostream &os_quality, std::ostream &os_cands, const my_queries &qs, const my_references &refs, const scoring_results &res, size_t pad, const bool ref_gaps, const papara_score_parameters &sp ) ;
    
    static void align_best_scores_oa( output_alignment *os, const my_queries &qs, const my_references &refs, const scoring_results &res, size_t pad, const bool ref_gaps, const papara_score_parameters &sp );
            
};



namespace {

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


template<typename state_t>
void gapstream_to_alignment( const std::vector<uint8_t> &gaps, const std::vector<state_t> &raw, std::vector<state_t> *out, state_t gap_char, const ref_gap_collector &rgc ) {

    typename std::vector<state_t>::const_reverse_iterator rit = raw.rbegin();


    size_t ref_ptr = rgc.ref_len();
    size_t gaps_left = rgc.gaps_before(ref_ptr);


    // this is kind of a hack: the 'inserted characters' are collected in insert, such that the additional
    // gaps can be put after the insert (remember, the trace is backward...). So the common inserts
    // are filled from left to right.
    // Except for the gap in the end (or more precisely in the beginning), where the gaps are put before the
    // insert, to make it appear as if the insert in the beginning is filled form right to left which looks
    // better.

    // the 'clean' solution would be to do 'de-novo' multiple alignment of the QS characters inside the common gaps...
    // TODO: do this and sell papara as a fragment assembler ;-)

    std::vector<state_t> insert;

    for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {

        if( *git == 1 || *git == 0 ) {
            // if a ref character is consumed, put in the collected inserts plus additional gaps if necessary.

            for( size_t i = 0; i < gaps_left; ++i ) {
                out->push_back(gap_char);
            }
            std::copy(insert.begin(), insert.end(), back_inserter(*out) );
            insert.clear();


            --ref_ptr;
            gaps_left = rgc.gaps_before(ref_ptr);
        }

        if ( *git == 1) {
            out->push_back(gap_char);

        } else if ( *git == 0 ) {
            assert( rit < raw.rend() );
            out->push_back(*rit);
            ++rit;

        } else {
            assert( rit < raw.rend() );
            //out->push_back(*rit);

            insert.push_back( *rit);
            ++rit;

            assert( gaps_left > 0 );
            --gaps_left;
        }
    }

    // NOTE: the insert comes before the gaps in this case
    std::copy(insert.begin(), insert.end(), back_inserter(*out) );
    for( size_t i = 0; i < gaps_left; ++i ) {
        out->push_back(gap_char);
    }

    std::reverse( out->begin(), out->end() );
}

template<typename state_t>
void gapstream_to_alignment_no_ref_gaps( const std::vector<uint8_t> &gaps, const std::vector<state_t> &raw, std::vector<state_t> *out, state_t gap_char ) {

    typename std::vector<state_t>::const_reverse_iterator rit = raw.rbegin();


    for ( std::vector<uint8_t>::const_iterator git = gaps.begin(); git != gaps.end(); ++git ) {

        if ( *git == 1) {
            out->push_back(gap_char);

        } else if ( *git == 0 ) {
            assert( rit < raw.rend() );
            out->push_back(*rit);
            ++rit;

        } else {
            assert( rit < raw.rend() );
            //out->push_back(*rit);

            ++rit;
        }
    }


    std::reverse( out->begin(), out->end() );
}


// void gapstream_to_position_map( const std::vector< uint8_t >& gaps, std::vector< int > &map) {
//     align_utils::trace_to_position_map( gaps, &map );
// 
// }

    
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
    return is.good();
}
} // end anonymous namespace for inline util functions (move to impl file!)



} // end of namespace papara

#endif

