#include <algorithm>
#include <deque>
#include <iomanip>

//#include <boost/bind.hpp>
#include "align_vec.h"
#define USE_BOOST_THREADS
#ifdef USE_BOOST_THREADS
#define BOOST_LIB_DIAGNOSTIC
#include <boost/thread.hpp>

namespace timpl = boost;

#else
#include "ivymike/thread.h"
namespace timpl = ivy_mike;
#endif

#include <boost/multi_array.hpp>
#include "ivymike/time.h"
#include "ivymike/write_png.h"

#include <pthread.h>

template <size_t W, typename seq_char_t>
struct db_block {
    int didx[W];
//     std::vector<seq_char_t> *ddata[W];
    size_t dpad[W];    
    size_t maxlen;
    int lj;
};

template <typename block>
struct block_queue {
    std::deque<block> m_blocks;
    timpl::mutex m_mtx;
    size_t m_ncup;
    block_queue() : m_ncup(0) {}
};

template <typename block_t>
struct worker {
    block_queue<block_t> &m_queue;
    
    worker( block_queue<block_t> &q ) : m_queue(q) {}
    void operator()() {
        
    }
};

template <size_t W, typename seq_char_t, typename score_t, typename sscore_t>
struct lworker {
    typedef db_block<W, seq_char_t> block_t;
    const int m_nthreads;
    const int m_rank;
    block_queue<block_t> &m_queue;
    const scoring_matrix &m_sm;
    const std::vector< std::vector<uint8_t> > &m_seq;
    const sscore_t gap_open;
    const sscore_t gap_extend;
    
    
    boost::multi_array<int,2> &m_outscore;
    
    lworker( int nthreads, int rank, block_queue<block_t>&q, const scoring_matrix &sm, const std::vector< std::vector<uint8_t> > &seq_, const sscore_t gap_open_, const sscore_t gap_extend_, boost::multi_array<int,2>&outscore ) 
    : m_nthreads(nthreads), m_rank(rank), m_queue(q), m_sm(sm), m_seq(seq_), gap_open(gap_open_), gap_extend(gap_extend_), m_outscore(outscore) {}
    
    void operator()() {
        
        size_t n_qseq = 0;
        size_t n_qchar = 0;
        
        size_t n_dseq = 0;
        size_t n_dchar = 0;
        
        bool first_block = true;
        aligned_buffer<seq_char_t> ddata_int;
        persistent_state<score_t> ps;
    
//         {
//             cpu_set_t cs;
//             CPU_ZERO( &cs );
//             CPU_SET( 0, &cs );
//             if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cs) != 0)
//             {
//                 printf("\n\nThere was a problem finding a physical core for thread number %d to run on.\n", 0);
//                 
//                 assert(0);
//             }
//             
//             
//         }
        
//         for( int i = 0; i < m_queue.m_blocks.size(); i++ ) {
//             
//             if( (i % m_nthreads) != m_rank ) {
//                 continue;
//             }
//         
//         
//             block_t block = m_queue.m_blocks[i];
        size_t ncups = 0;
        while(true) {
            block_t block;
            {
                timpl::lock_guard<timpl::mutex> lock( m_queue.m_mtx );
                
                if ( m_queue.m_blocks.empty() ) {
                    break;
                }
                block = m_queue.m_blocks.front();
                m_queue.m_blocks.pop_front();
            }
            
            
            //         std::cout << "sdis: " << sdi.size() << "\n";
//             std::cout << "asize: " << m_sm.num_states() << " " << block.maxlen << std::endl;
            aligned_buffer<sscore_t> qprofile( block.maxlen * W * m_sm.num_states());
            sscore_t *qpi = qprofile.begin();
            
            // setup the qprofile (= lookup table for match penalties along the db-sequences in the current block)
            // this is the faster (at least on core i5) two-step version, using interleaved db-sequences
            
            // setup buffer for interleaved db sequences
            if ( ddata_int.size() < block.maxlen * W ) {
                ddata_int.resize(block.maxlen * W);
            }
            
            // copy individual db sequences into interleaved buffer (padding the shorter sequnences
            seq_char_t *dint_iter = ddata_int.begin();
            const int zero_state = m_sm.get_zero_state();
            for ( int i = 0; i < block.maxlen; i++ ) {
                for ( int j = 0; j < W; j++ ) {
                    const std::vector<seq_char_t> &sdi = m_seq[block.didx[j]];//*(block.ddata[j]);
                    if ( i < sdi.size() ) {
                        *dint_iter = sdi[i];
                        
                    } else {
                        *dint_iter = zero_state;
                    }
                    
                    
                    
                    //                 std::cout << j << " " << int(*dint_iter) << " " << (i < sdi.size()) << "\n";
                    ++dint_iter;
                }
            }
            
            //copy interleaved scoring-matrix
            for ( int j = 0; j < m_sm.num_states(); j++ ) {
                dint_iter = ddata_int.begin();
                const char *cslice = m_sm.get_cslice(j);
                for ( int k = 0; k < block.maxlen; k++ ) {
                    for ( int l = 0; l < W; l++ ) {
                        //                     if( *dint_iter == zero_state ) {
                            //                         std::cout << int(cslice[*dint_iter]) << "\n";
                            //
                            //                     }
                            
                            *qpi = cslice[*dint_iter];
                            ++dint_iter;
                            ++qpi;
                    }
                }
            }
            
            std::vector<int> out(W);
            
            const size_t i_max = block.didx[block.lj];
//             const size_t i_max = m_seq.size() - 1;
            
            for ( size_t i_seq2 = 0; i_seq2 <= i_max; ++i_seq2 ) {
                const std::vector<uint8_t> &qdata = m_seq[i_seq2];
                
                if ( first_block ) {
                    n_qseq++;
                    n_qchar+=qdata.size();
                }
                
                
                
                align_vec<score_t,sscore_t,W>( ps, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                
                
                
                for ( int j = 0; j <= block.lj; j++ ) {
                    //                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
                    //                     std::cout << out[j] << "\t" << block.didx[j] << " " << i_seq2 << "\n";
//                     m_outscore[block.didx[j]][i_seq2] = out[j];
                    m_outscore[i_seq2][block.didx[j]] = out[j];
                    ncups += m_seq[i_seq2].size() * m_seq[block.didx[j]].size();
                }
                
                
                
            }
            
            for ( int j = 0; j <= block.lj; j++ ) {
                    //                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
//                     std::cout << out[j] << "\t" << block.didx[j] << " " << i_seq2 << "\n";
                    
                n_dseq++;
                n_dchar += m_seq[block.didx[j]].size();
            }
            first_block = false;
        }
        
        {
            std::cerr << n_qchar << " x " << n_dchar << "\n";
            timpl::lock_guard<timpl::mutex> lock( m_queue.m_mtx );
            m_queue.m_ncup += ncups;
        }
        
    }
};

void write_phylip_distmatrix( const boost::multi_array<int,2> &ma, const std::vector<std::string> &names, std::ostream &os ) {
    if( names.size() != ma.size() || ma.size() != ma[0].size() ) {
        throw std::runtime_error( "distance matrix seems fishy" );
    }
    os << ma.size() << "\n";
    os << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    for( int i = 0; i < ma.size(); i++ ) {
        os << names[i] << "\t";
        for( int j = 0; j < ma.size(); j++ ) {
            const float norm = std::max( ma[i][i], ma[j][j] );
            
            int mae;
            if( i <= j ) {
                mae = ma[i][j];
//                 mae = ma[j][i];
            } else {
                mae = ma[j][i];

            }
            
            const float dist = 1.0 - (mae / norm);
            
            os << dist << "\t";
        }
        os << "\n";
    }
}


void pairwise_seq_distance( const std::vector<std::string> &names, std::vector< std::vector<uint8_t> > &seq_raw, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread ) {
 #if 1
    const int W = 8;
    typedef short score_t;
    typedef short sscore_t;
#else
    const int W = 16;
    typedef unsigned char score_t;
    typedef char sscore_t;
#endif
    
//     size_t db_size = (sd.names.size() / W ) * W;
    

    ivy_mike::timer t1;
    
    
    std::vector< std::vector<uint8_t> > seq( seq_raw.size() );
//     seq.resize(400);
    for( int i = 0; i < seq.size(); i++ ) {
        std::for_each( seq_raw[i].begin(), seq_raw[i].end(), scoring_matrix::valid_state_appender<std::vector<uint8_t> >(sm, seq[i]) );
    }
    
//     const sscore_t gap_open = -5;
//     const sscore_t gap_extend = -2;
    
    typedef uint8_t seq_char_t;
    
    
   
    //std::string dname[W];
    
    
    
//     std::vector<score_t> dmask[W];
    
    
    
    bool have_input = true;
    
    
    size_t i_seq1 = 0;
    
    
    
    block_queue<db_block<W, seq_char_t> > q;
    std::deque<db_block<W, seq_char_t> > &blocks = q.m_blocks;
    while( have_input ) {
        
        
        
        
        
        // determine db sequences for the current block
        db_block<W, seq_char_t> block;
        block.maxlen = 0;
        
        block.lj = -1;
        for( int j = 0; j < W; j++ ) {
//             dname[j].resize(0);
            //ddata[j].resize(0);
//             have_input = (i_seq1 != seq.size());                       
           have_input = (i_seq1 != seq.size());
//             have_input = i_seq1 < 30;
           // std::cout << "have_input " << have_input << " " << seq.size() << "\n";
            
            
            
            // if there aren't enough db sequences left to fill the block, pad with last db sequence
            if( !have_input ) {
                
                // break immediately if there are no db sequences left (means #db-seqs % W == 0, or otherwise have_input would have been == false from last iteration)
                if( j == 0 ) {
                    break;
                } else {
//                     block.ddata[j] = block.ddata[block.lj];
                    block.didx[j] = block.didx[block.lj];
                }
            } else {
                block.didx[j] = i_seq1;
//                 block.ddata[j] = &seq[i_seq1];
                ++i_seq1;
                
                
                block.lj = j; // store largest valid 'j'
                
//                 for( int i = 0; i < ddata[j].length(); i++ ) {
//                         
//                     ddata[j][i] = sm.state_backmap(ddata[j][i]);
//                 }
            }
            
//             dmask[j].clear();
//             dmask[j].resize(ddata[j].length(), 0xffff );
            
            
            block.maxlen = std::max( block.maxlen, seq[block.didx[j]].size() );
        }
//         std::cout << "maxlen; " << block.maxlen << "\n";
        // jf == -1 at this point means that the block is empty (#db-seqs % W == 0)
        if( block.lj == -1 ) {
            break;
        }
       
        blocks.push_back(block);
    }
    
    std::cerr << "blocks: " << blocks.size() << "\n";
    //throw std::runtime_error( "exit" );
    
 
    
    
    
    boost::multi_array<int,2> out_scores;
    out_scores.resize( boost::extents[seq.size()][seq.size()] );
    
        
    
    
    timpl::thread_group tg;
    
    while( tg.size() < n_thread ) {
        lworker<W, seq_char_t, score_t, sscore_t> lw( n_thread, tg.size(), q, sm, seq, gap_open, gap_extend, out_scores );
        
        std::cerr << "thread " << tg.size() << "\n";
        
        tg.create_thread( lw );
    }
        
    tg.join_all();
    
    
    std::cerr << "aligned " << seq.size() << " x " << seq.size() << " sequences. " << q.m_ncup << " " << (q.m_ncup / (t1.elapsed() * 1.0e9)) << " GCup/s\n";
    
//     write_phylip_distmatrix( out_scores, names, std::cout );
    
    for( int i = 0; i < seq.size(); i++ ) {
        for( int j = 0; j < seq.size(); j++ ) {
            std::cout << out_scores[i][j] << "\t";
        }
        std::cout << "\n";
    }
    
//     ivy_mike::write_png( out_scores, std::cout );
}
