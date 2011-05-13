#include <fstream>
#include <memory>
#include <deque>
#include <boost/program_options.hpp>
#include "fasta.h"
#include "ivymike/write_png.h"
#include <ivymike/statistics.h>
#include "ivymike/thread.h"

void pairwise_seq_distance( const std::vector<std::string> &names, std::vector< std::vector<uint8_t> > &seq_raw, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );

class bla {
public:
    void operator()() {
        std::cout << "running\n";
    }
};

int main( int argc, char *argv[] ) {
    
//     bla b;
//     
//     //ivy_mike::thread t(b);
//     ivy_mike::thread_group tg;
//     
//     while( tg.size() < 2 ) {
//         bla b;
//         
//         tg.create_thread(b);
//     }
//     tg.join_all();
//     return 0;

    std::deque<float> v1;
    std::deque<float> v2;
    
    for( int i = 0; i < 1000; i++ ) {
        v1.push_back(sin(i));
        v2.push_back(sin(i + 3.1415/2));
    }
//     std::reverse( v2.begin(), v2.end() );
//     v2.push_back(4000);
    
    float cor = ivy_mike::correlation(v1, v2);
    std::cout << "cor: " << cor << "\n";
    
//    return 0;
    
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
        ("gap-extend,o", po::value<int>(&opt_gap_extend)->default_value(-3), "gap extend score (default: -3)")
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
    
    std::vector<std::string> qs_names;
    std::vector<std::vector<uint8_t> > qs_seqs;
    
    if( vm.count( "seq-file" ) != 1 ) {
        std::cout << desc << "\n";
        std::cout << "missing seq file\n";
        
        return -1;
    }
    
    std::ifstream qsf( opt_seq_file.c_str() );
    if( !qsf.good() ) {
        std::cout << "cannot open sequence file: " << opt_seq_file << "\n";
        
        return -1;
    }
    
    read_fasta( qsf, qs_names, qs_seqs);
    
    std::cerr << "using " << opt_threads << " threads\n";

    pairwise_seq_distance(qs_names, qs_seqs, *sm, opt_gap_open, opt_gap_extend, opt_threads);
    
}
