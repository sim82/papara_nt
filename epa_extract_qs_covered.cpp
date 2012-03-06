#include <iostream>

#include "tree_utils.h"
#include "ivymike/tree_parser.h"
#include "ivymike/multiple_alignment.h"
#include "ivymike/flat_map.h"
#include "boost/dynamic_bitset.hpp"

using ivy_mike::tree_parser_ms::ln_pool;
using ivy_mike::tree_parser_ms::parser;
using ivy_mike::tree_parser_ms::lnode;


bool is_gap( char c ) {
    return c == '-' || c == '?' || toupper(c) == 'N';
}

void pad( std::ostream &os, size_t n, char c ) {
    for( size_t i = 0; i < n; ++i ) {
        os << c;
    }
}

int main( int argc, char *argv[] ) {
    if( argc != 3 ) {
        std::cerr << "usage: " << argv[0] << "<tree file> <phylip file>\n";
        return 0;
    }
    
    ln_pool pool;
    
    parser p(argv[1], pool );
    lnode *n = p.parse();
    
    std::vector<std::string> tip_names;
    
    
    apply_lnode(n, [&](lnode *x) {
        if( x->m_data->isTip ) {
            tip_names.push_back( x->m_data->tipName );
        }
    });
    
    ivy_mike::multiple_alignment ma;
    ma.load_phylip(argv[2]);
    
    std::sort( tip_names.begin(), tip_names.end() );
    
    
    boost::dynamic_bitset<> bs_all(ma.data.front().size() );
    //bs_all.flip();
    
    size_t max_name_len = 0;
    for( size_t i = 0, e = ma.names.size(); i != e; ++i ) {
        const std::string &name = ma.names[i];
        max_name_len = std::max( max_name_len, name.size() );
        
        const auto &data = ma.data[i];
        
        if( !std::binary_search( tip_names.begin(), tip_names.end(), name ) ) {
            boost::dynamic_bitset<> bs(data.size());
            
            for( size_t j = 0; j < data.size(); ++j ) {
                bs[j] = !is_gap(data[j]);
            }
            
            
            bs_all |= bs;
        }
    }
    
    std::cout << ma.names.size() << " " << bs_all.count() << "\n";
    
    std::vector<size_t> unmasked;
    {
        size_t i = bs_all.find_first();
        while( i != bs_all.npos ) {
            unmasked.push_back(i);
            i = bs_all.find_next(i);
        }
    }
    for( size_t i = 0, e = ma.names.size(); i != e; ++i ) {
        std::cout << ma.names[i];
        const auto &data = ma.data[i];
        pad( std::cout, max_name_len - ma.names[i].size() + 1, ' ' );
        
        // excract 'unmasked' characters from current sequence
        std::transform(unmasked.begin(), unmasked.end(), std::ostream_iterator<char>(std::cout), [&](size_t u) { return data[u]; });
        std::cout << "\n";
    }
    
}