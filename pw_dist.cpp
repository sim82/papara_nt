#include <fstream>
#include <memory>
#include <deque>
#include <boost/program_options.hpp>
#include "fasta.h"
#include "ivymike/write_png.h"
#include <ivymike/statistics.h>
#include "ivymike/thread.h"
#include "ivymike/getopt.h"
#include <sys/mman.h>
#include <functional>
void pairwise_seq_distance( const std::vector<std::string> &names, std::vector< std::vector<uint8_t> > &seq_raw, const scoring_matrix &sm, const int gap_open, const int gap_extend, const int n_thread );

class bla {
public:
    void operator()() {
        std::cout << "running\n";
    }
};

int main( int argc, char *argv[] ) {

    
    

    ivy_mike::getopt::parser igp;
    
    igp.add_opt('h', false );
    igp.add_opt('t', true );
    igp.add_opt('f', true );
    igp.add_opt('m', true );
    igp.add_opt('n', true );
    igp.add_opt('o', true );
    igp.add_opt('e', true );
    igp.add_opt('s', true );
    igp.add_opt('t', true );
    
    
    igp.parse(argc, argv);
    
    if( igp.opt_count('h') != 0 ) {
        std::cout << 
        "  -h             print help message\n" <<
        "  -f arg         input sequence file (fasta)\n" <<
        "  -m arg         match score (implies DNA data, excludes option -s)\n" << 
        "  -n arg         mismatch score\n" << 
        "  -o arg         gap open score (default: -5, negtive means penalize)\n" <<
        "  -e arg         gap extend score (default: -3)\n" <<
        "  -s arg         scoring matrix (optional)\n" <<
        "  -t arg         number of threads (default: 1)\n\n" <<
        "The algorithm doesn't distinguish between DNA and AA data, as long as the\n" <<
        "input sequences are consistent with the scoring matrix. The use of the -m\n" <<
        "and -n options implies a DNA scoring matrix and will only work with DNA data\n";
        return 0;
        
    }
    
    
    if( igp.opt_count('f') != 1 ) {
        std::cerr << "missing option -f\n";
        return 0;
    }
    
    std::string opt_seq_file = igp.get_string('f');
    
    std::auto_ptr<scoring_matrix>sm;

    
    
    if( igp.opt_count('s') != 0 ) {
        std::string sm_name = igp.get_string('s');
        
        if( igp.opt_count('m') != 0 || igp.opt_count('n') != 0 ) {
            std::cerr << "option -s used in combination with option -m or -n\n";
        }
        
        
        std::ifstream is ( sm_name.c_str() );
        
        if( !is.good() ) {
            std::cout << "could not open scoring matrix " << sm_name << "\n";
            return -1;
        }
        
        sm.reset( new scoring_matrix( is ) );   
        
    } else {
        int match = 3;
        int mismatch = 0;
        
        igp.get_int_if_present('m',match);
        igp.get_int_if_present('n',mismatch);
        
        sm.reset( new scoring_matrix(match, mismatch));  
    }
    
//     std::cout << "file: " << igp.get_string('f');
    
    int opt_gap_open = -5;
    igp.get_int_if_present( 'o', opt_gap_open );
    int opt_gap_extend = -3;
    igp.get_int_if_present( 'e', opt_gap_extend );
  
    int opt_threads = 1;
    igp.get_int_if_present( 't', opt_threads );
    
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
    
    read_fasta( qsf, qs_names, qs_seqs);
    
    std::cerr << "using " << opt_threads << " threads\n";

    pairwise_seq_distance(qs_names, qs_seqs, *sm, opt_gap_open, opt_gap_extend, opt_threads);
    
}
