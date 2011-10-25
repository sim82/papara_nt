#include "ivymike/multiple_alignment.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <map>
#include <functional>
#include <cstring>
#include <boost/io/ios_state.hpp>

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/bind.hpp>
#include <boost/array.hpp>



#include "parsimony.h"
#include "pvec.h"

#include "pars_align_seq.h"
#include "pars_align_gapp_seq.h"
#include "fasta.h"
#include "vec_unit.h"
#include "align_pvec_vec.h"

#include "ivymike/tree_parser.h"
#include "ivymike/time.h"
#include "ivymike/getopt.h"
#include "ivymike/thread.h"
#include "ivymike/demangle.h"
#include "ivymike/stupid_ptr.h"
#include "ivymike/algorithm.h"
#include "ivymike/smart_ptr.h"


using namespace ivy_mike::tree_parser_ms;


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




//class my_adata : public ivy_mike::tree_parser_ms::adata {
////     static int ct;
//    //std::vector<parsimony_state> m_pvec;
//    pvec_t m_pvec;
//
//public:
////     int m_ct;
//    my_adata_gen() {
//
////         std::cout << "my_adata\n";
//
//    }
//
//    virtual ~my_adata_gen() {
//
////         std::cout << "~my_adata\n";
//
//    }
//
//    virtual void visit() {
////         std::cout << "tr: " << m_ct << "\n";
//    }
//    void init_pvec(const std::vector< uint8_t >& seq) {
//
//
//        m_pvec.init( seq );
////         std::cout << "init_pvec: " << m_pvec.size() << "\n";
////                 m_pvec.reserve(seq.size());
////         for( std::vector< uint8_t >::const_iterator it = seq.begin(); it != seq.end(); ++it ) {
////             m_pvec.push_back(dna_to_parsimony_state(*it));
////
////         }
//    }
//    pvec_t &get_pvec() {
//        return m_pvec;
//    }
//
//
//};


class my_ldata : public ldata {

public:
	typedef std::vector<double> pvec_t;



private:
	boost::array<pvec_t,4> pvecs_;
};

class my_fact : public ivy_mike::tree_parser_ms::node_data_factory {

    virtual my_ldata *alloc_ldata() {

        return new my_ldata;
    }

};

class references {
public:
	references( sptr::shared_ptr<ln_pool> pool, const std::string &tree_name, const std::string &ali_name )
	  : pool_(pool), tree_name_(tree_name), ali_name_(ali_name)
	{}


	void preprocess() {


	}

private:
	sptr::shared_ptr<ln_pool> pool_;

	const std::string tree_name_;
	const std::string ali_name_;

};


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
    bool opt_use_gpu;
    bool opt_force_overwrite;

    igp.add_opt( 't', igo::value<std::string>(opt_tree_name) );
    igp.add_opt( 's', igo::value<std::string>(opt_alignment_name) );
    igp.add_opt( 'q', igo::value<std::string>(opt_qs_name) );
    igp.add_opt( 'c', igo::value<bool>(opt_use_cgap, true).set_default(false) );
    igp.add_opt( 'j', igo::value<int>(opt_num_threads).set_default(1) );
    igp.add_opt( 'n', igo::value<std::string>(opt_run_name).set_default("default") );
    igp.add_opt( 'b', igo::value<bool>(opt_write_testbench, true).set_default(false) );
    igp.add_opt( 'g', igo::value<bool>(opt_use_gpu, true).set_default(false) );
    igp.add_opt( 'f', igo::value<bool>(opt_force_overwrite, true).set_default(false) );

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

    if( opt_run_name != "default" && !opt_force_overwrite && file_exists(log_filename.c_str()) ) {
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

    lout << "bla\n";
    //return 0;


    sptr::shared_ptr<ln_pool> pool(new ln_pool(std::auto_ptr<node_data_factory>(new my_fact<my_adata>)) );


    references refs( pool, opt_tree_name, opt_alignment_name );

    if( !opt_use_cgap ) {
        pnt_ptr.reset( new papara_nt<pvec_pgap, 4>( opt_tree_name.c_str(), opt_alignment_name.c_str(), qs_name, opt_write_testbench, opt_use_gpu ));
    } else {
        pnt_ptr.reset( new papara_nt<pvec_cgap, 8>( opt_tree_name.c_str(), opt_alignment_name.c_str(), qs_name, opt_write_testbench, opt_use_gpu ));
    }

//     return 0;

    std::cerr << "using " << opt_num_threads << " threads\n";


    papara_nt_i &pnt = *pnt_ptr;
    pnt.calc_scores( opt_num_threads );

    {
        std::ofstream os( filename( opt_run_name, "scores" ).c_str() );
        pnt.print_best_scores(os);
    }

    {
        std::ofstream os( filename( opt_run_name, "alignment" ).c_str() );
        std::ofstream os_quality( filename( opt_run_name, "quality" ).c_str() );

        //         pnt.dump_ref_seqs(os);
        //         pnt.align_best_scores(os);
        pnt.write_result_phylip(os, os_quality);
    }

    {
    	std::ofstream os( filename( opt_run_name, "all_scores" ).c_str() );
    	pnt.write_all_scores(os);
    }

    //ivymike::LN *n = tp.parse();

//     getchar();
    //ivymike::LN::free( n );
//     delete n;


    std::cout << t.elapsed() << std::endl;
    lout << "SUCCESS " << t.elapsed() << std::endl;
    return 0;
//     getchar();
}

