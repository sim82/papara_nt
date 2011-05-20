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


#include <fstream>
#include <memory>
#include <deque>

#include "fasta.h"
#include "ivymike/write_png.h"
#include <ivymike/statistics.h>
#include "ivymike/thread.h"
#include "ivymike/getopt.h"

#include <functional>
#include <iomanip>
// #define PWDIST_INLINE
#include "pairwise_seq_distance.h"
//void pairwise_seq_distance( const std::vector<std::string> &names, std::vector< std::vector<uint8_t> > &seq_raw, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );





void write_phylip_distmatrix( const ivy_mike::tdmatrix<int> &ma, const std::vector<std::string> &names, std::ostream &os ) {
    if( names.size() != ma.size() || ma.size() != ma[0].size() ) {
        throw std::runtime_error( "distance matrix seems fishy" );
    }
    os << ma.size() << "\n";
    os << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    for( int i = 0; i < ma.size(); i++ ) {
        os << names[i] << "\t";
        for( int j = 0; j < ma.size(); j++ ) {
            
            // three modes for normalizing: min, max and mean
            //const float norm = std::min( ma[i][i], ma[j][j] );
//             const float norm = std::max( ma[i][i], ma[j][j] );
            const float norm = (ma[i][i] + ma[j][j]) * 0.5;
            
            
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


class bla {
public:
    void operator()() {
        std::cout << "running\n";
    }
};

int main( int argc, char *argv[] ) {

//     tdmatrix<int> tdm( 10, 10 );
// 
//     std::fill( tdm.begin(), tdm.end(), 33 );
//     
//     tdm[5][1] = 10;
//     
//     for( tdmatrix<int>::row_iterator rit = tdm.row_begin(); rit != tdm.row_end(); ++rit ) {
//         odmatrix<int> odm = *rit;
//         for( int i = 0; i < odm.size(); i++ ) {
//             std::cout << odm[i] << " ";
//         }
//         std::cout << "\n";
//     }
//     
//     return 0;
    
    ivy_mike::getopt::parser igp;
    std::string opt_seq_file;
    int opt_match;
    int opt_mismatch;
    int opt_gap_open;
    int opt_gap_extend;
    std::string opt_sm_name;
    int opt_threads;
    
    bool opt_out_dist_matrix;
    bool opt_out_score_matrix;
    bool opt_out_pgm_image;
    bool opt_out_none;
    
    igp.add_opt('h', false );
    
    igp.add_opt('f', ivy_mike::getopt::value<std::string>(opt_seq_file) );
    igp.add_opt('m', ivy_mike::getopt::value<int>(opt_match).set_default(3) );
    igp.add_opt('n', ivy_mike::getopt::value<int>(opt_mismatch).set_default(0) );
    igp.add_opt('o', ivy_mike::getopt::value<int>(opt_gap_open).set_default(-5) );
    igp.add_opt('e', ivy_mike::getopt::value<int>(opt_gap_extend).set_default(-3) );
    igp.add_opt('s', ivy_mike::getopt::value<std::string>(opt_sm_name) );
    igp.add_opt('t', ivy_mike::getopt::value<int>(opt_threads).set_default(1) );
    igp.add_opt('1', ivy_mike::getopt::value<bool>(opt_out_dist_matrix, true).set_default(false) );
    igp.add_opt('2', ivy_mike::getopt::value<bool>(opt_out_score_matrix, true).set_default(false) );
    igp.add_opt('3', ivy_mike::getopt::value<bool>(opt_out_pgm_image, true).set_default(false) );
    igp.add_opt('4', ivy_mike::getopt::value<bool>(opt_out_none, true).set_default(false) );
    
    bool ret = igp.parse(argc, argv);
    
    
    if( !opt_out_dist_matrix && !opt_out_score_matrix && !opt_out_pgm_image && !opt_out_none) {
        opt_out_dist_matrix = true;
    }
    
//     std::cout << "opt_match: " << &opt_seq_file << "\n";
    
//     return 0;
    
    if( igp.opt_count('h') != 0 || !ret ) {
        std::cout << 
        "  -h        print help message\n" <<
        "  -f arg    input sequence file (fasta)\n" <<
        "  -m arg    match score (implies DNA data, excludes option -s)\n" << 
        "  -n arg    mismatch score\n" << 
        "  -o arg    gap open score (default: -5, negtive means penalize)\n" <<
        "  -e arg    gap extend score (default: -3)\n" <<
        "  -s arg    scoring matrix (optional)\n" <<
        "  -t arg    number of threads (default: 1)\n\n" <<
        "  -1        output distance matrix (PHYLIP format, e.g. for nj-tree building with ninja)\n" <<
        "  -2        output raw score matrix\n" <<
        "  -3        output greyscale pgm image (gimmick)\n" <<
        "  -4        output no results (e.g., for benchmark)\n" <<
        " In any case, the output will be written to stdout.\n\n" <<
        "The algorithm doesn't distinguish between DNA and AA data, as long as the\n" <<
        "input sequences are consistent with the scoring matrix. The use of the -m\n" <<
        "and -n options implies a DNA scoring matrix and will only work with DNA data\n";
        return 0;
        
    }
    
    
    if( igp.opt_count('f') != 1 ) {
        std::cerr << "missing option -f\n";
        return 0;
    }
    
   // std::string opt_seq_file = igp.get_string('f');
    
    std::auto_ptr<scoring_matrix>sm;

    
    
    if( igp.opt_count('s') != 0 ) {
        //std::string sm_name = igp.get_string('s');
        
        if( igp.opt_count('m') != 0 || igp.opt_count('n') != 0 ) {
            std::cerr << "option -s used in combination with option -m or -n\n";
        }
        
        
        std::ifstream is ( opt_sm_name.c_str() );
        
        if( !is.good() ) {
            std::cout << "could not open scoring matrix " << opt_sm_name << "\n";
            return -1;
        }
        
        std::cerr << "using generic scoring matrix from file: " << opt_sm_name << "\n";
        
        sm.reset( new scoring_matrix( is ) );   
        
    } else {
//         int match = 3;
//         int mismatch = 0;
//         
//         igp.get_int_if_present('m',match);
//         igp.get_int_if_present('n',mismatch);
        
        std::cerr << "using flat DNA scoring. match: " << opt_match << " mismatch: " << opt_mismatch << "\n";
        
        sm.reset( new scoring_matrix(opt_match, opt_mismatch));  
    }
    
//     std::cout << "file: " << igp.get_string('f');
    
    std::cerr << "gap open  : " << opt_gap_open << "\n";
    std::cerr << "gap extend: " << opt_gap_extend << "\n";
    std::cerr << "nthreads  : " << opt_threads << "\n";
    std::cerr << "seq. file : " << opt_seq_file << "\n";
    
#if 0    
//     return 0;
     
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    
    std::string opt_seq_file;
    int opt_match;
    int opt_mismatch;
    int opt_gap_open;
    int opt_gap_extend;
    int opt_threads;
    std::string opt_scoring_matrix;
    
    desc.add_options()
        ("help,h", "print help message" )
        ("seq-file,f", po::value<std::string>(&opt_seq_file), "input sequence file (fasta)" )
        ("match,m", po::value<int>(&opt_match)->default_value(3), "match score" )
        ("mismatch,n", po::value<int>(&opt_mismatch)->default_value(0), "mismatch score" )
        ("gap-open,o", po::value<int>(&opt_gap_open)->default_value(-5), "gap open score (default: -5, negtive means penalize)")
        ("gap-extend,e", po::value<int>(&opt_gap_extend)->default_value(-3), "gap extend score (default: -3)")
        ("scoring-matrix,s", po::value<std::string>(&opt_scoring_matrix), "scoring matrix (optional)")
        ("threads,t", po::value<int>(&opt_threads)->default_value(1), "number of threads (default: 1)" );
//     po::positional_options_description p;
//     p.add( "sdf-file", -1 );
    
    
    po::variables_map vm;
   // po::store( po::command_line_parser( argc, argv ), desc.positional(p).run(), vm );
    try {
        po::store( po::parse_command_line( argc, argv, desc ), vm );
    } catch( po::error x ) {
     
        std::cout << "could not parse commanline: " << x.what() << "\n";
        std::cout << "available options:\n" << desc << "\n";
        return -1;
        
    }
    
    po::notify(vm);
    
    if( vm.count("help" ) ) {
        std::cout << desc << "\n";
        return -1;
    }
    
    
    std::auto_ptr<scoring_matrix>sm;

    if( vm.count( "scoring-matrix" ) == 1 ) {
        if( vm.count( "match" ) != 0 || vm.count( "mismatch" ) != 0 ) {
            std::cout << "command line error: give either explicite scores OR scoring matrix\n";
            std::cout << desc << "\n";
            
            return -1;
        }
        
        std::ifstream is ( opt_scoring_matrix.c_str() );
        
        if( !is.good() ) {
            std::cout << "could not open scoring matrix " << opt_scoring_matrix << "\n";
            return -1;
        }
        
        sm.reset( new scoring_matrix( is ) );   
    } else {
        sm.reset( new scoring_matrix(opt_match, opt_mismatch));
    }
    if( vm.count( "seq-file" ) != 1 ) {
        std::cout << desc << "\n";
        std::cout << "missing seq file\n";
        
        return -1;
    }
    
#endif
    std::vector<std::string> qs_names;
    std::vector<std::vector<uint8_t> > qs_seqs;
    
    
    
    std::ifstream qsf( opt_seq_file.c_str() );
    if( !qsf.good() ) {
        std::cout << "cannot open sequence file: " << opt_seq_file << "\n";
        
        return -1;
    }
    
    read_fasta( qsf, *sm, qs_names, qs_seqs);
    
    std::cerr << "using " << opt_threads << " threads\n";
//     return 0;
    
    //     write_phylip_distmatrix( out_scores, names, std::cout );
    ivy_mike::tdmatrix<int> out_scores( qs_seqs.size(), qs_seqs.size() );
    pairwise_seq_distance( qs_seqs, out_scores, *sm, opt_gap_open, opt_gap_extend, opt_threads);
    
    if( opt_out_dist_matrix ) {
        write_phylip_distmatrix( out_scores, qs_names, std::cout );
    } else if( opt_out_score_matrix ) {
        for( int i = 0; i < qs_seqs.size(); i++ ) {
            for( int j = 0; j < qs_seqs.size(); j++ ) {
                std::cout << out_scores[i][j] << "\t";
            }
            std::cout << "\n";
        }
    } else if( opt_out_pgm_image ) {
        ivy_mike::write_png( out_scores, std::cout );        
    }

    
}
