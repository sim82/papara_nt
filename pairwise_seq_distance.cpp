#include <algorithm>
#include <deque>
#include <boost/bind.hpp>
#include "align_vec.h"
#include <boost/thread.hpp>
#include "ivymike/time.h"


template <size_t W, typename seq_char_t>
struct db_block {
    int didx[W];
    std::vector<seq_char_t> *ddata[W];
    size_t dpad[W];    
    size_t maxlen;
    int lj;
};

template <typename block>
struct block_queue {
    std::deque<block> m_blocks;
    boost::mutex m_mtx;
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
    block_queue<block_t> &m_queue;
    const scoring_matrix &m_sm;
    const std::vector< std::vector<uint8_t> > &m_seq;
    const sscore_t gap_open;
    const sscore_t gap_extend;
    const int m_nthreads;
    const int m_rank;
    
    lworker( int nthreads, int rank, block_queue<block_t>&q, scoring_matrix &sm, const std::vector< std::vector<uint8_t> > &seq_, const sscore_t gap_open_, const sscore_t gap_extend_ ) 
    : m_nthreads(nthreads), m_rank(rank), m_queue(q), m_sm(sm), m_seq(seq_), gap_open(gap_open_), gap_extend(gap_extend_) {}
    
    void operator()() {
        
        size_t n_qseq = 0;
        size_t n_qchar = 0;
        
        size_t n_dseq = 0;
        size_t n_dchar = 0;
        
        bool first_block = true;
        aligned_buffer<seq_char_t> ddata_int;
        persistent_state<score_t> ps;
    
        
        
//         for( int i = 0; i < m_queue.m_blocks.size(); i++ ) {
//             
//             if( (i % m_nthreads) != m_rank ) {
//                 continue;
//             }
//         
//         
//             block_t block = m_queue.m_blocks[i];
        
        while(true) {
            block_t block;
            {
                boost::lock_guard<boost::mutex> lock( m_queue.m_mtx );
                
                if ( m_queue.m_blocks.empty() ) {
                    break;
                }
                block = m_queue.m_blocks.front();
                m_queue.m_blocks.pop_front();
            }
            
            
            //         std::cout << "sdis: " << sdi.size() << "\n";
            std::cout << "asize: " << m_sm.num_states() << " " << block.maxlen << std::endl;
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
                    std::vector<seq_char_t> &sdi = *(block.ddata[j]);
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
            
            
            
            for ( size_t i_seq2 = 0; i_seq2 < m_seq.size(); ++i_seq2 ) {
                const std::vector<uint8_t> &qdata = m_seq[i_seq2];
                
                if ( first_block ) {
                    n_qseq++;
                    n_qchar+=qdata.size();
                }
                
                
                
                align_vec<score_t,sscore_t,W>( ps, block.maxlen, qdata, m_sm, qprofile, gap_open, gap_extend, out );
                
                for ( int j = 0; j <= block.lj; j++ ) {
                    //                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
//                     std::cout << out[j] << "\t" << block.didx[j] << " " << i_seq2 << "\n";
                    
             
                }
                
                
                
            }
            
            for ( int j = 0; j <= block.lj; j++ ) {
                    //                 std::cout << out[j] << "\t" << dname[j] << " " << qname << " " << ddata[j].size() << "\n";
//                     std::cout << out[j] << "\t" << block.didx[j] << " " << i_seq2 << "\n";
                    
                n_dseq++;
                n_dchar += block.ddata[j]->size();
            }
            first_block = false;
        }
        
        {
            std::cout << n_qchar << " x " << n_dchar << "\n";
            boost::lock_guard<boost::mutex> lock( m_queue.m_mtx );
            m_queue.m_ncup += n_qchar * n_dchar;
        }
        
    }
};

void pairwise_seq_distance( std::vector< std::vector<uint8_t> > &seq_raw ) {
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
    
    scoring_matrix sm( 3, 0 );
    std::vector< std::vector<uint8_t> > seq( seq_raw.size() );
//     seq.resize(400);
    for( int i = 0; i < seq.size(); i++ ) {
        std::for_each( seq_raw[i].begin(), seq_raw[i].end(), scoring_matrix::valid_state_appender<std::vector<uint8_t> >(sm, seq[i]) );
    }
    
    const sscore_t gap_open = -5;
    const sscore_t gap_extend = -2;
    
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
            
            have_input = (i_seq1 != seq.size());
//             have_input = i_seq1 < 30;
            std::cout << "have_input " << have_input << " " << seq.size() << "\n";
            
            
            
            // if there aren't enough db sequences left to fill the block, pad with last db sequence
            if( !have_input ) {
                
                // break immediately if there are no db sequences left (means #db-seqs % W == 0, or otherwise have_input would have been == false from last iteration)
                if( j == 0 ) {
                    break;
                } else {
                    block.ddata[j] = block.ddata[block.lj];
                    block.didx[j] = -1;
                }
            } else {
                block.didx[j] = i_seq1;
                block.ddata[j] = &seq[i_seq1];
                ++i_seq1;
                
                
                block.lj = j; // store largest valid 'j'
                
//                 for( int i = 0; i < ddata[j].length(); i++ ) {
//                         
//                     ddata[j][i] = sm.state_backmap(ddata[j][i]);
//                 }
            }
            
//             dmask[j].clear();
//             dmask[j].resize(ddata[j].length(), 0xffff );
            
            
            block.maxlen = std::max( block.maxlen, block.ddata[j]->size() );
        }
        
        // jf == -1 at this point means that the block is empty (#db-seqs % W == 0)
        if( block.lj == -1 ) {
            break;
        }
        blocks.push_back(block);
    }
    //throw std::runtime_error( "exit" );
    
 
    
    
    
    
    
        
    const size_t n_thread = 4;
    
    boost::thread_group tg;
    
    while( tg.size() < n_thread ) {
        lworker<W, seq_char_t, score_t, sscore_t> lw( n_thread, tg.size(), q, sm, seq, gap_open, gap_extend );
        
        tg.create_thread( lw );
    }
        
    tg.join_all();
    
    
    std::cout << "aligned " << seq.size() << " x " << seq.size() << " sequences. " << q.m_ncup << " " << (q.m_ncup / (t1.elapsed() * 1.0e9)) << " GCup/s\n";
    
}