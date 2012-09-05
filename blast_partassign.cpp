#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <cctype>
#include <cassert>
#include <map>

#include "blast_partassign.h"

using partassign::partition;
using partassign::blast_hit;


blast_hit partassign::next_hit( std::istream &is ) {
    std::string line;
    
    assert( is.good() );
    
    std::getline( is, line );
    
    std::istringstream ss(line);
    // TODO: ok returning an default constructed blast_hit to signal EOF is not the best idea... I actually planned that EOF 
    // should be checked _before_ calling this function, but it didn't turn out to work...
    
    if( line.empty() ) {
        return blast_hit();
    }
    
//     for( size_t i = 0; i < 12; ++i ) {
//         int s;
//         ss >> s;
//         std::cout << i << " " << s << "\n";
//         
//         
//     }
    
    std::string qs, ref;
    float ident;
    int len, mismatch, gap_open, qs_start, qs_end, ref_start, ref_end;
    std::string evalue;
    blast_hit hit;
    
    bool valid = ss >> hit.qs_name >> hit.ref_name >> ident >> len >> mismatch >> gap_open >> hit.qs_start >> hit.qs_end >> hit.ref_start >> hit.ref_end >> evalue >> hit.bit_score;
    
    // convert indices into proper form
    hit.ref_start -= 1;
    hit.ref_end -= 1;
    hit.qs_start -= 1;
    hit.qs_end -= 1;
    
    if( !valid ) {
        std::cerr << "bad line: " << line << "\n";
        throw std::runtime_error( "could not read line in blast output file.\n" );
    }
    
    return hit;
    
}

template<typename T>
T from_string( const std::string &str ) {
    std::istringstream ss(str);
    T x;
    ss >> x;
    return x;
}
// static bool not_space( char x ) {
//     return !isspace(x); // doing the same with pre-c++11 functional is ridiculus
// }

partition partassign::next_partition( std::istream &is ) {
    // TODO: there is absolutely no (=ZILCH) error checking in this function
    
    std::string line;
    
    assert( is.good() );
    
    std::getline( is, line );
    
    std::istringstream ss(line);

    // TODO: d.t.o.
    if( line.empty() ) {
        return partition();
    }
    
    
    // parse model name:
    // it = start if line (maybe scan for non-space?)
    // it_next = first ','
    std::string::iterator it = line.begin();
    std::string::iterator it_next = std::find( it, line.end(), ',' );
    
    std::string model_name( it, it_next );
//     std::cout << "model: '" << model_name << "'\n";
    
    it = it_next + 1;
    
    
//     std::find_if( it, line.end(), not_space );
    
    // parse gene name:
    // it = first non_space character after the ','
    // it_next = first space after the gene name
    it = std::find_if( it, line.end(), std::not1(std::ptr_fun<int,int>(std::isspace) ) );
    it_next = std::find_if( it, line.end(), std::ptr_fun<int,int>(std::isspace) );
    
    std::string gene_name( it, it_next );
//     std::cout << "gene: '" << gene_name << "'\n";
    
    it = std::find( it_next, line.end(), '=' );
    ++it;
    it = std::find_if( it, line.end(), std::not1(std::ptr_fun<int,int>(std::isspace) ) );
    
    it_next = std::find( it, line.end(), '-' );
    
    std::string start_str( it, it_next );
    it = it_next + 1;
    std::string end_str( it, line.end() );
//     std::cout << "start: '" << start_str << "'\n";
//     std::cout << "end: '" << end_str << "'\n";
//     
    partition part;
    part.start = from_string<int>(start_str) - 1;
    part.end = from_string<int>(end_str) - 1;
    return part;
}

// int main() {
//     
//     {
//         std::ifstream is( "test.model" );
//         assert( is.good() );
//         
//         while( is.good() ) {
//             partition p = partassign::next_partition(is);
//             
//             if( p.start == -1 ) {
//                 break;
//             }
//             
//             std::cout << "part: " << p.start << " " << p.end << "\n";
//             
//             
//         }
//         
//     }
//     
//     {
//         std::ifstream is( "test.out" );
//         assert( is.good() );
//         
//         std::map<std::string,blast_hit> hit_map;
//         
//         while( is.good() ) {
//             
//             blast_hit hit = partassign::next_hit( is );
//             
//             if( hit.qs_name.empty() ) {
//                 break;
//             }
//             
//             //std::cout << "hit: " << hit.qs_start << " " << hit.bit_score << "\n";
//             
//             
//             // check for multiple hits per QS, keep the hit with the highest bit-score
//             std::map< std::string, blast_hit >::iterator it = hit_map.lower_bound( hit.qs_name );
//             if( it != hit_map.end() && it->second.qs_name == hit.qs_name ) {
//                 
//                 if( it->second.bit_score < hit.bit_score ) {
//                     it->second = hit;
//                     //                 std::cout << "replace\n";
//                 }
//             } else {
//                 hit_map.insert( it, std::make_pair( hit.qs_name, hit ));
//             }
//         }
//     }
//     
// }
namespace partassign
{
part_assignment::part_assignment ( std::istream& blast_out, std::istream& part_file )
{
    while ( part_file.good() ) {
        partitions_.push_back ( partassign::next_partition ( part_file ) );
        if ( partitions_.back().start == -1 ) { // returning a partition with negaitve indices is next_partition's way of signalling EOF
            partitions_.pop_back();
            break;
        }
    }


    std::map<std::string,partassign::blast_hit> hit_map;

    while ( blast_out.good() ) {

        blast_hit hit = partassign::next_hit ( blast_out );

        if ( hit.qs_name.empty() ) {
            break;
        }

        //std::cout << "hit: " << hit.qs_start << " " << hit.bit_score << "\n";


        // check for multiple hits per QS, keep the hit with the highest bit-score
        std::map< std::string, blast_hit >::iterator it = hit_map.lower_bound ( hit.qs_name );
        if ( it != hit_map.end() && it->second.qs_name == hit.qs_name ) {

            if ( it->second.bit_score < hit.bit_score ) {
                it->second = hit;
                //                 std::cout << "replace\n";
            }
        } else {
            hit_map.insert ( it, std::make_pair ( hit.qs_name, hit ) );
        }
    }
    hits_.swap(hit_map);

//     for ( std::map<std::string,blast_hit>::iterator it = hit_map.begin(), eit = hit_map.end(); it != eit; ++it ) {
//         int start = it->second.ref_start;
//         int end = it->second.ref_end;
// 
//         int part_idx = -1;
// 
//         for ( size_t i = 0; i < partitions_.size(); ++i ) {
//             const partassign::partition &part = partitions_[i];
// 
//             if ( start >= part.start && end <= part.end ) {
//                 part_idx = int ( i );
//                 break;
//             }
//         }
// 
//         if ( part_idx == -1 ) {
//             std::cerr << "QS cannot be uniqely assigned to a single partition: " << it->first << "[" << start << "-" << end << "\n";
//             throw std::runtime_error ( "partitons incompatible with balst hits" );
//         }

//         a [it->first] = part_idx;
// 
//     }

}
const blast_hit& part_assignment::get_blast_hit ( const std::string& qs_name ) const 
{
    std::map< std::string, partassign::blast_hit >::const_iterator it = hits_.find ( qs_name );
    if ( it == hits_.end() ) {
        throw std::runtime_error ( "qs name not found" );
    }

    return it->second;
}
// const partition& part_assignment::partition ( const std::string& name ) const
// {
//     std::map< std::string, int >::const_iterator it = assignments_.find ( name );
// 
//     if ( it == assignments_.end() ) {
//         std::stringstream ss;
//         ss << "no partition assignment for qs " << name;
//         throw std::runtime_error ( ss.str() );
//     }
// 
//     return partitions_.at ( it->second );
// }
}
