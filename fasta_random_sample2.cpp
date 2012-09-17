/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */

// FIXME: someone tries to sabotage us on windows by sneaking in the awful min/max macros. Strangely it only happens in this source file... 



#include "ivymike/disable_shit.h" 
#include <iostream>
#include <algorithm>
#include <iterator>



#include "ivymike/fasta.h"
#include "stepwise_align.h"

void process_fasta( std::istream &is )
{

	
    std::vector<std::string> primers;
    primers.push_back("CTTGGTCATTTAGAGGAAGT");
    primers.push_back("CGATGAAGAACGCAG");
    
    std::vector<std::string> names;
    std::vector<std::vector<uint8_t> > data;
    
    ivy_mike::read_fasta( is, names, data, false );
    
    const size_t n = names.size();
    assert( n == data.size());
    
    std::vector<size_t> rind;
    for( size_t i = 0; i < names.size(); ++i ) {
        rind.push_back(i);
    }
    
    std::random_shuffle( rind.begin(), rind.end() );
    
    const size_t target_size = 10;
    
    size_t num_written = 0;
    
    
    int match_score = 5;
    ivy_mike::scoring_matrix sm(5, -3);
    
    
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
		return -1;
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

