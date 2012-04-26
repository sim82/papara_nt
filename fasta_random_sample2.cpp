#include <iostream>
#include <algorithm>
#include <iterator>

#include "fasta.h"
#include "stepwise_align.h"


void process_fasta( std::istream &is )
{
    std::vector<std::string> primers;
    primers.push_back("CTTGGTCATTTAGAGGAAGT");
    primers.push_back("CGATGAAGAACGCAG");
    
    std::vector<std::string> names;
    std::vector<std::vector<uint8_t> > data;
    
    read_fasta( is, names, data, false );
    
    const size_t n = names.size();
    assert( n == data.size());
    
    std::vector<size_t> rind;
    for( size_t i = 0; i < names.size(); ++i ) {
        rind.push_back(i);
    }
    
    std::random_shuffle( rind.begin(), rind.end() );
    
    const size_t target_size = 5;
    
    size_t num_written = 0;
    
    
    int match_score = 5;
    scoring_matrix sm(5, -3);
    
    
    while( num_written < target_size && !rind.empty()) {
        size_t r = rind.back();
        rind.pop_back();
        
        
        
        std::vector<uint8_t> seq = data.at(r);
        
        
        
        size_t num_rejects = 0;
        for( size_t j = 0; j < primers.size(); ++j ) {
            std::vector<uint8_t> p( primers[j].begin(), primers[j].end() );
            float score = align_freeshift( sm, seq, p, -5, -2, false );
            float escore = match_score * p.size();
            if( score < 0.9 * escore ) {
                ++num_rejects;
            }
            
//             std::cout << "score " << j << ": " << score << " " << p.size() * match_score << "\n";
        }
        
        if( num_rejects != 0 ) {
//             std::cout << "reject!\n";
            continue;
        }
        
        std::cout << ">" << names.at(r) << "\n";
        
        std::copy( data.at(r).begin(), data.at(r).end(), std::ostream_iterator<char>(std::cout) );
        std::cout << "\n";        
        
        ++num_written;
    }
    
    
}

int main( int argc, char *argv[] ) {
 
    
    if( argc < 2 ) {
        std::cerr << "missing random seed/inout files" << std::endl;
    }
    
    unsigned int seed = atoi( argv[1] );
    std::cerr << "rand seed: " << seed << std::endl;
    srand(seed);
    
    if( argc == 2 ) {
        process_fasta( std::cin );
    } else {
        for( int i = 2; i < argc; ++i ) {
            std::ifstream is(argv[i] );
            
            if( !is.good() ) {
                std::cerr << "WARNING: cannot read: " << argv[i] << "\n";
                continue;
            }
            
            process_fasta( is );
        }
    }
}

