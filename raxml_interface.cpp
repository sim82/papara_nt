#include <Poco/Process.h>
#include <Poco/File.h>
#include <Poco/Pipe.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <deque>
#include <vector>

#include "ivymike/time.h"

#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tr1/unordered_map.hpp>
#include <boost/tr1/array.hpp>
#include "raxml_interface.h"
#include "tree_utils.h"
#include "tree_similarity.h"


using ivy_mike::tree_parser_ms::lnode;

namespace {


boost::shared_ptr<ivy_mike::tree_parser_ms::ln_pool> g_pool;

}


static uint8_t transform_to_gap( uint8_t c ) {
	if( c == '-' || std::toupper(c) == 'N' ) {
		return '1';
	} else {
		return '0';
	}
}

static void write_phylip( const std::map<std::string,const std::vector<uint8_t> * const> &name_to_seq, std::ostream &os, bool gap_model ) {
	size_t len = size_t(-1);
	size_t max_name_len = 0;


	for( std::map<std::string,const std::vector<uint8_t> * const>::const_iterator it = name_to_seq.begin(), e = name_to_seq.end(); it != e; ++it ) {

		if( len == size_t(-1) ) {
			len = it->second->size();
		} else {
			if( len != it->second->size() ) {
				throw std::runtime_error( "non-uniform sequence lengths in name_to_seq map" );
			}
		}

		max_name_len = std::max( max_name_len, it->first.size() );
	}


	os << name_to_seq.size() << " " << len << "\n";
	for( std::map<std::string,const std::vector<uint8_t> * const>::const_iterator it = name_to_seq.begin(), e = name_to_seq.end(); it != e; ++it ) {
		os << std::setw(max_name_len + 1) << std::left << it->first;
		if( !gap_model ) {
			std::copy( it->second->begin(), it->second->end(), std::ostream_iterator<uint8_t>(os));
		} else {
			std::transform( it->second->begin(), it->second->end(), std::ostream_iterator<uint8_t>(os), &transform_to_gap );
		}
		os << "\n";
	}



}

void pipe_into_deque( Poco::Pipe &pipe, std::deque<char> &deq ) {
	// the deque version might be much more efficient if _lots_ of data are read from the pipe...

	// i don't trust the complicated implementation for std::vector, so use this...

	//std::vector<char> buf(4096);
	std::tr1::array<char,4096> buf;
	size_t readsize = -1;

	while( (readsize = pipe.readBytes( buf.data(), buf.size())) != 0 ) {
		deq.insert( deq.end(), buf.begin(), buf.begin() + readsize );
	}
}

//void pipe_into_vector( Poco::Pipe &pipe, std::vector<char> &vec ) {
//	const size_t incsize = 4096;
//	vec.resize(incsize);
//	std::vector<char>::iterator start_it = vec.begin();
//	std::vector<char>::iterator end_it = vec.end();
//	size_t readsize = -1;
//	while( (readsize = pipe.readBytes(&(*start_it), end_it - start_it)) != 0 ) {
//		//std::copy( raxbuf.begin(), raxbuf.begin() + readsize, std::ostream_iterator<char>(std::cout) );
//		start_it += readsize;
//		assert( start_it <= end_it );
//
//
//
//		if( start_it == end_it ) {
//			// increase buffer size
//			const size_t oldsize = vec.size();
//			vec.resize(oldsize + incsize );
//
//
//			// recreate start/end iterators
//			start_it = vec.begin() + oldsize;
//			end_it = vec.end();
//
////			std::cout << "resize: " << raxbuf.size() << "\n";
//		}
//
//		assert( start_it < end_it );
//	}
//
//	vec.resize( start_it - vec.begin() );
//}

void optimize_branch_lengths( ivy_mike::tree_parser_ms::lnode *tree, const std::map<std::string,const std::vector<uint8_t> * const> &name_to_seq ) {
	ivy_mike::perf_timer perf_timer(!true);

	std::vector< std::pair< ivy_mike::tree_parser_ms::lnode*, ivy_mike::tree_parser_ms::lnode* > > edges;
	std::vector<boost::dynamic_bitset<> > splits;
	std::vector<ivy_mike::tree_parser_ms::lnode *> sorted_tips;



	get_all_splits( tree, edges, splits, sorted_tips );

//	for( std::vector<boost::dynamic_bitset<> >::iterator it = splits.begin(); it != splits.end(); ++it ) {
//		std::cout << "split: " << (*it) << "\n";
//	}
//
//	for( std::vector<lnode*>::iterator it = sorted_tips.begin(); it != sorted_tips.end(); ++it ) {
//		//(*it)->m_data->print(std::cout);
//		//std::cout << "\n";
//
//		std::cout << "tip: " << (*(*it)->m_data) << "\n";
//
//	}


	assert( edges.size() == splits.size() );

	std::cout << "edges: " << edges.size() << "\n";

	if( edges.size() < 4 ) {

		return;
	}


	perf_timer.add_int();
	const bool gap_model = false;
	{
		std::ofstream ost( "tmp_tree" );
		ivy_mike::tree_parser_ms::print_newick( next_non_tip( towards_tree( tree )), ost );


		std::ofstream osa( "tmp_ali" );
		write_phylip( name_to_seq, osa, gap_model );
	}


	perf_timer.add_int();
	Poco::Process::Args args;
	args.push_back( "-T" );
	args.push_back( "4" );
	args.push_back( "-f" );
	args.push_back( "e" );
	args.push_back( "-m" );
	if( gap_model ) {
		args.push_back( "BINCAT" );
	} else {
		args.push_back( "GTRCAT" );
	}
	args.push_back( "-n" );
	args.push_back( "Txxx" );
	args.push_back( "-s" );
	args.push_back( "tmp_ali" );
	args.push_back( "-t" );
	args.push_back( "tmp_tree" );
	args.push_back( "-n" );
	args.push_back( "Txxx" );

	std::string raxml( "/home/sim/src_exelixis/standard-RAxML/raxmlHPC-PTHREADS-SSE3" );


	// remove old raxml output
	Poco::File f( "RAxML_info.Txxx" );

	if( f.exists() ) {
		assert( f.exists() && f.canWrite() && f.isFile() && "file not readable" );
		f.remove();
	}


#if 1
	Poco::Pipe raxout_pipe;

	Poco::ProcessHandle proc = Poco::Process::launch( raxml, args, 0, &raxout_pipe, 0 );


//	std::vector<char> raxbuf;
//	pipe_into_vector( raxout_pipe, raxbuf );
	std::deque<char> raxbuf;
	pipe_into_deque( raxout_pipe, raxbuf );

	std::cout << "raxml wrote " << raxbuf.size() << "\n";
#else
	Poco::ProcessHandle proc = Poco::Process::launch( raxml, args, 0, 0, 0 );
#endif

	int ret = proc.wait();

	perf_timer.add_int();
	std::cout << "wait for raxml: " << ret << "\n";


	{
		if( g_pool.get() == 0 ) {
			g_pool.reset( new ivy_mike::tree_parser_ms::ln_pool() );
		}
		const char *raxml_tree = "RAxML_result.Txxx";
		ivy_mike::tree_parser_ms::parser p( raxml_tree, *g_pool );

		ivy_mike::tree_parser_ms::lnode *rax_tree = p.parse();

		tip_collector_dumb<ivy_mike::tree_parser_ms::lnode> tc;
		visit_lnode(rax_tree, tc );

		std::cout << "rax tree: " << tc.m_nodes.size() << "\n";

		std::vector< std::pair<lnode*,lnode*> > rax_edges;
		std::vector<boost::dynamic_bitset<> > rax_splits;
		std::vector<ivy_mike::tree_parser_ms::lnode *> rax_sorted_tips;

		get_all_splits( rax_tree, rax_edges, rax_splits, rax_sorted_tips );

//		for( std::vector<lnode*>::iterator it = rax_sorted_tips.begin(); it != rax_sorted_tips.end(); ++it ) {
//			//(*it)->m_data->print(std::cout);
//			//std::cout << "\n";
//
//			std::cout << "raxtip: " << (*(*it)->m_data) << "\n";
//
//		}


//		for( std::vector<boost::dynamic_bitset<> >::iterator it = rax_splits.begin(); it != rax_splits.end(); ++it ) {
//			std::cout << "rax split: " << (*it) << "\n";
//		}


		std::vector<std::pair<boost::dynamic_bitset<>, std::pair<lnode*,lnode*> > > sorted_splits;


		// hmm hmm, maybe tr1::unordered_map is the better alternative... test!
		{
			assert( rax_edges.size() == edges.size() );
			assert( rax_edges.size() == rax_splits.size() );

			sorted_splits.resize( rax_edges.size() );
			for( size_t i = 0; i < rax_edges.size(); ++i ) {
				sorted_splits[i].first.swap(rax_splits[i]);
				sorted_splits[i].second = rax_edges[i];
			}

			std::sort( sorted_splits.begin(), sorted_splits.end() );
			for( size_t i = 0; i < sorted_splits.size(); ++i ) {
				rax_splits[i].swap(sorted_splits[i].first);
				rax_edges[i] = sorted_splits[i].second;
			}
		}


		for( size_t i = 0; i < splits.size(); ++i ) {


//			std::vector<boost::dynamic_bitset<> >::iterator it = std::lower_bound(rax_splits.begin(), rax_splits.end(), splits[i] );

			std::vector<boost::dynamic_bitset<> >::iterator it = std::find( rax_splits.begin(), rax_splits.end(), splits[i] );
//			std::cout << "find: " << itf - rax_splits.begin() << " " << it - rax_splits.begin() << "\n";

		//	assert( it != )

			if( it == rax_splits.end() ) {
				throw std::runtime_error( "split not found in original tree" );
			} else if( (*it) != splits[i] ) {
				throw std::runtime_error( "split not found in original tree (2)" );
			}

			size_t rax_idx = it - rax_splits.begin();
			std::pair<lnode*,lnode*> edge = edges[i];

			assert( rax_idx < rax_edges.size() && "lookup in rax_splits produced out of range index for rax_edges");
			std::pair<lnode*,lnode*> rax_edge = rax_edges[i];

			assert( rax_edge.first->backLen == rax_edge.second->backLen && "inconsitstent backLen values in rax tree" );
			edge.first->backLen = rax_edge.first->backLen;
			edge.second->backLen = rax_edge.second->backLen;


		}

		//for( size_t i = 0; < )


		g_pool->clear();
		g_pool->sweep();

	}

	perf_timer.add_int();

	perf_timer.print();
}
