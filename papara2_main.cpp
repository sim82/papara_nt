#include "blast_partassign.h"
#include "ivymike/getopt.h"
#include "ivymike/time.h"

#include "papara.h"

using namespace papara;

namespace papara {
log_stream lout;
}

void print_banner( std::ostream &os ) {
    os << "______    ______    ______        ___\n";
    os << "| ___ \\   | ___ \\   | ___ \\   |\\ | |   \n";
    os << "| |_/ /_ _| |_/ /_ _| |_/ /__ | \\| | \n";
    os << "|  __/ _` |  __/ _` |    // _` |\n";
    os << "| | | (_| | | | (_| | |\\ \\ (_| |\n";
    os << "\\_|  \\__,_\\_|  \\__,_\\_| \\_\\__,_|\n";
    os << "  Version 2.2\n";

}

void pad( std::ostream &os, size_t w ) {
    for( size_t i = 0; i < w; ++i ) {
        os << " ";
    }
}

void print_help( std::ostream &os, const std::vector<std::string> &options, const std::vector<std::string> &text ) {
    assert( options.size() == text.size() );

    size_t max_opt_size = 0;

    for( size_t i = 0; i < options.size(); ++i ) {
        max_opt_size = std::max( max_opt_size, options[i].size() );
    }


    size_t pad_right = 3;
    size_t col_width = 2 + max_opt_size + pad_right;
    for( size_t i = 0; i < options.size(); ++i ) {
        pad( os, 2 );
        os << options[i];
        pad( os, col_width - 2 - options[i].size() );
        //os << text[i];
        const std::string &t = text[i];


        for( size_t j = 0; j < t.size(); ++j ) {
            if( t[j] == '@') {
                os << "\n";
                pad( os, col_width );
            } else {
                os << t[j];
            }
        }

        os << "\n";
    }
    
    
}


void print_help( std::ostream &os ) {
    std::vector<std::string> options;
    std::vector<std::string> text;

    options.push_back( "-t <tree file>" );
    text.push_back( "File name of the reference tree (newick format)");

    options.push_back( "-s <ref alignment>" );
    text.push_back( "File name of the reference alignment (phylip format)");

    options.push_back( "-q <query seqs.>" );
    text.push_back( "File name of the query sequences (fasta format)");

    options.push_back( "-a" );
    text.push_back( "Sequences are protein data");

    options.push_back( "-j <num threads>" );
    text.push_back( "Specify number of threads (default:1)");

    options.push_back( "-n <run name>" );
    text.push_back( "Specify filename suffix of the output@files (default: \"default\")");

    options.push_back( "-r" );
    text.push_back( "Turn of writing RA-side gaps in the output file.");

    options.push_back( "-p" );
    text.push_back( "User defined scoring scheme: <open>:<extend>:<match>:<match cg>@The default scores correspond to '-p -3:-1:2:-3'" );
    
    print_help( os, options, text );

}


template<typename pvec_t, typename seq_tag>
std::vector<std::pair<size_t,size_t> > resolve_qs_bounds( references<pvec_t,seq_tag> &refs, queries<seq_tag> &qs, const partassign::part_assignment &part_assign ) {
    std::vector<std::pair<size_t,size_t> > bounds;
    
    for( size_t i = 0; i < qs.size(); ++i ) {
        const std::string &qs_name = qs.name_at(i);
        const partassign::blast_hit &hit = part_assign.get_blast_hit( qs_name );
        
         size_t ref_idx = refs.find_name( hit.ref_name );
         
         if( ref_idx == size_t(-1) ) {
             throw std::runtime_error( "ref name of blast hit not found" );
         }
         const std::vector<int> &ng_map = refs.ng_map_at(ref_idx);
         
         if( size_t(hit.ref_start) >= ng_map.size() || size_t(hit.ref_end) >= ng_map.size() ) {
             throw std::runtime_error( "blast hit region outside of reference sequence" );
         }
         
         // map position in (non-gappy) ref sequence onto alignment column
         int col_start = ng_map.at(hit.ref_start);
         int col_end = ng_map.at(hit.ref_end);
         
         int part_idx = -1;
         
         std::vector< partassign::partition > partitions = part_assign.partitions();
         for ( size_t i = 0; i < partitions.size(); ++i ) {
             const partassign::partition &part = partitions[i];
             
             if ( col_start >= part.start && col_end <= part.end ) {
                 part_idx = int ( i );
                 break;
             }
         }
         
         if ( part_idx == -1 ) {
             std::cerr << "QS cannot be uniquely assigned to a single partition: " << qs_name << " [" << col_start << "-" << col_end << "]\n";
             throw std::runtime_error ( "partitons incompatible with blast hits" );
         }
         
         bounds.push_back( std::make_pair( partitions[part_idx].start, partitions[part_idx].end ));
         
    }

    
    return bounds;
}

template<typename pvec_t, typename seq_tag>
void run_papara( const std::string &qs_name, const std::string &alignment_name, const std::string &tree_name, size_t num_threads, const std::string &run_name, const bool ref_gaps, const papara_score_parameters &sp, bool write_fasta, partassign::part_assignment *part_assign ) {

    ivy_mike::perf_timer t1;

    queries<seq_tag> qs(qs_name.c_str());

    
    
    
    t1.add_int();
    references<pvec_t,seq_tag> refs( tree_name.c_str(), alignment_name.c_str(), &qs );

    
    
    t1.add_int();

    qs.preprocess();
   
    t1.add_int();

    refs.remove_full_gaps();
    refs.build_ref_vecs();

    if( part_assign != 0 ) {
        //qs.init_partition_assignments( *part_assign );
        std::vector<std::pair<size_t,size_t> > qs_bounds = resolve_qs_bounds( refs, qs, *part_assign );
        
        
        qs.set_per_qs_bounds( qs_bounds );
        
        
    }
    
    t1.add_int();

//     t1.print();

    const size_t num_candidates = 0;

    scoring_results res( qs.size(), scoring_results::candidates(num_candidates) );


    lout << "scoring scheme: " << sp.gap_open << " " << sp.gap_extend << " " << sp.match << " " << sp.match_cgap << "\n";

    driver<pvec_t,seq_tag>::calc_scores(num_threads, refs, qs, &res, sp );

    std::string score_file(filename(run_name, "alignment"));
    std::string quality_file(filename(run_name, "quality"));
    std::string cands_file(filename(run_name, "cands"));


    size_t pad = 1 + std::max(qs.max_name_length(), refs.max_name_length());

//     std::ofstream os( score_file.c_str() );
//     assert( os.good() );

    std::ofstream os_qual( quality_file.c_str() );
    assert( os_qual.good() );


    std::ofstream os_cands;
    if( num_candidates > 1 ) {
        os_cands.open( cands_file.c_str() );
        assert( os_cands.good());
    }



    std::auto_ptr<papara::output_alignment> oa;
    if( write_fasta ) {
        oa.reset( new papara::output_alignment_fasta( score_file.c_str() ));
    } else {
        oa.reset( new papara::output_alignment_phylip( score_file.c_str() ));
    }
    
    //refs.write_seqs(os, pad);
    //     driver<pvec_t,seq_tag>::align_best_scores( os, os_qual, os_cands, qs, refs, res, pad, ref_gaps, sp );
    driver<pvec_t,seq_tag>::align_best_scores_oa( oa.get(), qs, refs, res, pad, ref_gaps, sp );
    
}



int main( int argc, char *argv[] ) {

//     aligned_buffer<int> xxx(1024);
    

    
    namespace igo = ivy_mike::getopt;

    ivy_mike::getopt::parser igp;

    std::string opt_tree_name;
    std::string opt_alignment_name;
    std::string opt_qs_name;
    std::string opt_user_parameters;
    
    std::string opt_blast_hits;
    std::string opt_partitions;
    
    bool opt_use_cgap;
    int opt_num_threads;
    std::string opt_run_name;
    bool opt_write_testbench;
    bool opt_force_overwrite;
    bool opt_aa;
    bool opt_no_ref_gaps;
    bool opt_print_help;
    bool opt_write_fasta;
    
    igp.add_opt( 't', igo::value<std::string>(opt_tree_name) );
    igp.add_opt( 's', igo::value<std::string>(opt_alignment_name) );
    igp.add_opt( 'q', igo::value<std::string>(opt_qs_name) );
    igp.add_opt( 'c', igo::value<bool>(opt_use_cgap, true).set_default(false) );
    igp.add_opt( 'a', igo::value<bool>(opt_aa, true).set_default(false) );
    igp.add_opt( 'j', igo::value<int>(opt_num_threads).set_default(1) );
    igp.add_opt( 'n', igo::value<std::string>(opt_run_name).set_default("default") );
    igp.add_opt( 'b', igo::value<bool>(opt_write_testbench, true).set_default(false) );
    igp.add_opt( 'f', igo::value<bool>(opt_force_overwrite, true).set_default(false) );
    igp.add_opt( 'r', igo::value<bool>(opt_no_ref_gaps, true).set_default(false) );
    igp.add_opt( 'p', igo::value<std::string>(opt_user_parameters).set_default("") );
    igp.add_opt( 'h', igo::value<bool>(opt_print_help, true).set_default(false) );
    igp.add_opt( 'g', igo::value<bool>(opt_write_fasta, true).set_default(false) );
    igp.add_opt( 'l', igo::value<std::string>(opt_blast_hits) );
    igp.add_opt( 'x', igo::value<std::string>(opt_partitions) );
    
    igp.parse(argc,argv);

    if( opt_print_help ) {
        print_banner( std::cerr );
        std::cerr << "\nUser options:\n";
        print_help( std::cerr);
        return 0;
        
    }
         
    
    if( igp.opt_count('t') != 1 || igp.opt_count('s') != 1 ) {
        print_banner(std::cerr);

        std::cerr << "missing options -t and/or -s (-q is optional)\n";

        print_help( std::cerr );
        return 0;
    }
    
    // optional accelration by blast hits/partition file
    std::auto_ptr<partassign::part_assignment> part_assignment;
    
    if( igp.opt_count('l') == 1 || igp.opt_count('x') == 1 ) {
        if( igp.opt_count('l') != igp.opt_count('x') ) {
            std::cerr << "options -b and -x have to be used together (or not at all)\n";
            print_help( std::cerr );
            return 0;
        }
        
        std::ifstream blast_is( opt_blast_hits.c_str() );
        if( !blast_is.good() ) {
            std::cerr << "can not open blast hits file\n";
            return 0;
        }
        
        std::ifstream part_is( opt_partitions.c_str() );
        if( !part_is.good() ) {
            std::cerr << "can not open partition file\n";
            return 0;
        }
        
        
        part_assignment.reset( new partassign::part_assignment( blast_is, part_is ));
    }
    
    
    ivy_mike::timer t;

//    const char *qs_name = 0;
//    if( !opt_qs_name.empty() ) {
//        qs_name = opt_qs_name.c_str();
//    }
//
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
    
    

    const bool ref_gaps = !opt_no_ref_gaps;


    papara_score_parameters sp = papara_score_parameters::default_scores();
    if( !opt_user_parameters.empty() ) {
        std::cout << "meeeeep\n";
        
        sp = papara_score_parameters::parse_scores(opt_user_parameters.c_str());
    }
    
    
    if( opt_use_cgap ) {

        if( opt_aa ) {
            run_papara<pvec_cgap, tag_aa>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get() );
        } else {
            run_papara<pvec_cgap, tag_dna>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get() );
        }
    } else {
        if( opt_aa ) {
            run_papara<pvec_pgap, tag_aa>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get() );
        } else {
            run_papara<pvec_pgap, tag_dna>( opt_qs_name, opt_alignment_name, opt_tree_name, opt_num_threads, opt_run_name, ref_gaps, sp, opt_write_fasta, part_assignment.get() );
        }
    }

    std::cout << t.elapsed() << std::endl;
    lout << "SUCCESS " << t.elapsed() << std::endl;

    
    
    
    return 0;
//     getchar();
}

