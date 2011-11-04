#include <Poco/Process.h>
#include <Poco/File.h>
#include <Poco/Pipe.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <deque>
#include <vector>
#include <cassert>

#include "ivymike/time.h"

#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/tr1/unordered_map.hpp>
#include <boost/tr1/array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lexical_cast.hpp>
#include "raxml_interface.h"
#include "tree_utils.h"
#include "tree_similarity.h"


using ivy_mike::tree_parser_ms::lnode;
using ivy_mike::tree_parser_ms::adata;
using ivy_mike::tree_parser_ms::ln_pool;

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


lnode *optimize_branch_lengths2( ivy_mike::tree_parser_ms::lnode *tree, const std::map<std::string, const std::vector<uint8_t> * const> &name_to_seq, ln_pool &pool ) {
	ivy_mike::perf_timer perf_timer(!true);


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

	tip_collector_dumb<ivy_mike::tree_parser_ms::lnode> tc;
	visit_lnode(tree, tc );


	std::map<std::string, sptr::shared_ptr<lnode> > tip_map;

	for( std::vector<lnode *>::iterator it = tc.m_nodes.begin(); it != tc.m_nodes.end(); ++it ) {
		tip_map.insert( std::make_pair( (*it)->m_data->tipName, (*it)->get_smart_ptr() ) );
	}


	std::cout << "nodes: " << tc.m_nodes.size() << "\n";

	if( tc.m_nodes.size() < 4 ) {

		return 0;
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
	args.push_back( "Tyyy" );
	args.push_back( "-s" );
	args.push_back( "tmp_ali" );
	args.push_back( "-t" );
	args.push_back( "tmp_tree" );

	std::string raxml( "/home/sim/src_exelixis/standard-RAxML/raxmlHPC-PTHREADS-SSE3" );


	// remove old raxml output
	Poco::File f( "RAxML_info.Tyyy" );

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


	const char *raxml_tree = "RAxML_result.Tyyy";
	ivy_mike::tree_parser_ms::parser p( raxml_tree, pool );

	ivy_mike::tree_parser_ms::lnode *rax_tree = p.parse();

	tip_collector_dumb<ivy_mike::tree_parser_ms::lnode> rax_tc;
	visit_lnode(rax_tree, rax_tc );

	lnode *ret_node = rax_tree;
	// link adata from old tree into corresponding tips of the new tree.
	for( std::vector<lnode *>::iterator it = rax_tc.m_nodes.begin(); it != rax_tc.m_nodes.end(); ++it ) {


		std::string name = (*it)->m_data->tipName;

		std::map<std::string, sptr::shared_ptr<lnode> >::iterator op = tip_map.find( name );

		assert( op != tip_map.end() && "could not find tip from new tree in old tree");
//		sptr::shared_ptr<adata> old_adata( op->second->m_data );

		lnode *old_node = op->second.get();

		lnode *new_node = *it;

		if( new_node == ret_node ) {
			ret_node = old_node;
		}

		assert( old_node->back != 0 );
		assert( new_node->back != 0 );

		new_node->back->back = old_node;
		//new_node->back->backLen = old_node->backLen;

		old_node->back->back = 0;
		old_node->back = new_node->back;
		old_node->backLen = new_node->backLen;

		new_node->back = 0;
		//sptr::shared_ptr<lnode> old_node =
	}

	perf_timer.add_int();

	perf_timer.print();

	assert( ret_node->back != 0 );
	assert( ret_node->back->back != 0 );
	return ret_node;
}

namespace ublas = boost::numeric::ublas;

lnode *generate_marginal_ancestral_state_pvecs( ln_pool &pool, const std::string &tree_name, const std::string &ali_name, std::vector<ublas::matrix<double> > *pvecs ) {
	ivy_mike::perf_timer perf_timer(!true);



	perf_timer.add_int();


	perf_timer.add_int();
	Poco::Process::Args args;

	args.push_back( "-f" );
	args.push_back( "A" );
	args.push_back( "-m" );

	args.push_back( "GTRGAMMA" );

	args.push_back( "-n" );
	args.push_back( "Y1" );
	args.push_back( "-s" );
	args.push_back( ali_name );
	args.push_back( "-t" );
	args.push_back( tree_name );

	std::string raxml( "/home/sim/src_exelixis/hacked_as_raxml/raxmlHPC" );


	// remove old raxml output
	Poco::File f( "RAxML_info.Y1" );

	if( f.exists() ) {
		assert( f.exists() && f.canWrite() && f.isFile() && "file not readable" );
		f.remove();
	}


#if 0
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
#else
	int ret = 42;
#endif

	perf_timer.add_int();
	std::cout << "wait for raxml: " << ret << "\n";

	//assert( ret == 42 );


	const char *raxml_tree = "RAxML_nodeLabelledRootedTree.Y1";
	ivy_mike::tree_parser_ms::parser p( raxml_tree, pool );

	ivy_mike::tree_parser_ms::lnode *rax_tree = p.parse();

	std::ifstream pis( "RAxML_marginalAncestralProbabilities.Y1", std::ios::binary );
	assert( pis.good() );


	// read the binary ancestral probability file.
	// this is a huge improvement over the text file stuff (below) and now reads at 1200Mb/s
	// from warm disk-cache instead of 20 Mb/s...

	pvecs->clear();
	ivy_mike::timer t1;

	std::vector<double> tmp;
	std::string line;
	std::vector<double> token;
	std::stringstream strstr;

	while( !pis.eof() ) {
		int32_t counter;
		pis.read((char*)&counter, 4 );



		if( counter == -1 ) {
			break;
		}

		assert( size_t(counter) == pvecs->size() );

		int32_t width;
		pis.read((char*)&width, 4 );

//		std::cout << "width: " << width << "\n";
		assert( width > 0 );



		pvecs->push_back(ublas::matrix<double>(width, 4));
		ublas::matrix<double> &mat = pvecs->back();

		// the underlying unbounded_array seems to be guaranteed to have a sensible memory layout. read straight into it.
		// TODO: get fancy and directly back the boost::matrix with the mmaped binary file ;-)
		pis.read((char*)mat.data().begin(), width * 4 * 8 );
	}

	size_t s = -1;
	for( size_t i = 0; i < pvecs->size(); ++i ) {
		if( s == size_t(-1) ) {
			s = (*pvecs)[i].size1();
		}

		if( (*pvecs)[i].size1() != s ) {
			throw std::runtime_error( "inconsistent ancestral pvec lengths");
		}

	}

	std::cout << "read anc probs: " << t1.elapsed() << "\n";

	return rax_tree;

}


#if 0

// text file version. slow as hell
lnode *generate_marginal_ancestral_state_pvecs( ln_pool &pool, const std::string &tree_name, const std::string &ali_name, std::vector<ublas::matrix<double> > *pvecs ) {
	ivy_mike::perf_timer perf_timer(!true);



	perf_timer.add_int();


	perf_timer.add_int();
	Poco::Process::Args args;

	args.push_back( "-f" );
	args.push_back( "A" );
	args.push_back( "-m" );

	args.push_back( "GTRGAMMA" );

	args.push_back( "-n" );
	args.push_back( "X1" );
	args.push_back( "-s" );
	args.push_back( ali_name );
	args.push_back( "-t" );
	args.push_back( tree_name );

	std::string raxml( "/home/sim/src_exelixis/hacked_as_raxml/raxmlHPC" );


	// remove old raxml output
	Poco::File f( "RAxML_info.X1" );

	if( f.exists() ) {
		assert( f.exists() && f.canWrite() && f.isFile() && "file not readable" );
		f.remove();
	}


#if 0
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
#else
	int ret = 42;
#endif

	perf_timer.add_int();
	std::cout << "wait for raxml: " << ret << "\n";

	assert( ret == 42 );


	const char *raxml_tree = "RAxML_nodeLabelledRootedTree.X1";
	ivy_mike::tree_parser_ms::parser p( raxml_tree, pool );

	ivy_mike::tree_parser_ms::lnode *rax_tree = p.parse();

	std::ifstream pis( "RAxML_marginalAncestralProbabilities.X1" );
	assert( pis.good() );


	pvecs->clear();

	// while this loop is a bit ugly, it seems to be the fastest way to
	// read the state prob. text files. It still only parses about 20 Mb/s, because
	// string->double conversion is incredibly slow...

	// keeping the allocation of this stuff outside the loop helps a lot.
	std::vector<double> tmp;
	std::string line;
	std::vector<double> token;
	std::stringstream strstr;

	while( !pis.eof() ) {
		// read line and tokenize into doubles

		line.clear();
		std::getline( pis, line );

//		std::stringstream strstr(line);
		strstr.clear();
		strstr.str(line);

		std::istream_iterator<double> it(strstr);
		const static std::istream_iterator<double> end;

		token.assign(it, end);

		// depending on the number of token it is either a header/data/end-line
		if( token.size() == 1 ) {
			assert( pvecs->size() == size_t(token[0]) );
			tmp.clear();
		} else if( token.size() == 4 ) {
			std::copy( token.begin(), token.end(), std::back_inserter(tmp) );
		} else {
			assert( tmp.size() % 4 == 0 );
			const size_t nlines = tmp.size() / 4;

			pvecs->push_back( ublas::matrix<double>() );

			ublas::matrix<double> &mat = pvecs->back();
			mat.resize( nlines, 4 );

			// Using the copy should also work, as mat is a (nlines x 4) row-major matrix.
			// TODO: look if ublas makes any guarantees about the layout of matrix::data().
			//std::copy( tmp.begin(), tmp.end(), mat.data().begin() );

			std::vector<double>::iterator tit = tmp.begin();
			ublas::matrix<double>::iterator1 mit1 = mat.begin1();

			// fold the linear (row-major) tmp vector into the ublas::matrix.
			// FIXME: isn't there some way to wrap a std::vector in a ublas::matrix_expression?
			for( size_t i = 0; i < nlines; ++i, tit+=4, ++mit1 ) {
				assert( tit <= tmp.end() - 4 );
				assert( mit1 < mat.end1() );
				std::copy( tit, tit + 4, mit1.begin() );
			}
		}
	}

	size_t s = -1;
	for( size_t i = 0; i < pvecs->size(); ++i ) {
		if( s == size_t(-1) ) {
			s = (*pvecs)[i].size2();
		}

		if( (*pvecs)[i].size2() != s ) {
			throw std::runtime_error( "inconsistent ancestral pvec lengths");
		}

	}


	return rax_tree;

}
#endif
