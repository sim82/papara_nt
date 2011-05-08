/*
 * Copyright (C) 2010 Simon A. Berger
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

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <stdint.h>
#include <cstdlib>
#include <cstddef>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "aligned_buffer.h"




class mapped_file {
    int m_fd;
    
    void *m_base;
    size_t m_size;
    size_t m_ptr;
    
    mapped_file() {}
    mapped_file( const mapped_file &mf );
    const mapped_file &operator=(const mapped_file &other );
public:
    
    mapped_file( const char *name ) 
        : m_ptr(0)
    {
        m_fd = open( name ,O_RDONLY );
        
        if( m_fd == -1 ) {
         
            throw new std::runtime_error( "cannot open file for mapping" );
        }
        
        m_size = lseek( m_fd, 0, SEEK_END );
        
        
        m_base = mmap( 0, m_size, PROT_READ, MAP_SHARED, m_fd, 0 );
        
        if( m_base == 0 ) {
            throw std::runtime_error( "mmap failed" );   
        }
        
        madvise( m_base, m_size, MADV_SEQUENTIAL );
        
    }
    
    ~mapped_file() {
        munmap( m_base, m_size );
    }


    inline void check_bounds() {
     
        
    }
    int get() {
        if( eof() ) {
            return -1;
        } else {
            return ((char*)m_base)[m_ptr++];
        }
    }
 
    void unget() {
        if( m_ptr > 0 ) {
            m_ptr--;
        }
    }
 
    bool eof() {
        return m_ptr >= m_size;   
    }
    bool good() {
        return !eof();
    }
    
    void seekg( off_t ptr ) {
        m_ptr = ptr;   
    }
    void clear() {
        
    }
};


template<class input>
static inline void get_line( input &is, std::vector<char> &linebuf ) {
    linebuf.resize(0);
    while(true) {
        char c = is.get();
        if( c == '\n' || is.eof() ) {
            break;
        }
        linebuf.push_back(c);
    }

}

static inline bool xisspace( int c ) {
    return c == ' ' || c == '\n';   
}

template<class input>
static void read_fasta( input &is, std::vector<std::string> &names, std::vector<std::string> &data ) {
 
//     std::vector<char> linebuf;//1024 * 1024); // uhm, this is ugly...
    //std::string linebuf;
    
    names.clear();
    data.clear();
    std::string *data_accum = 0;
    
    while( is.good() && !is.eof() ) {
    
        //is.getline( linebuf.data(), linebuf.size() );
            
            
        int c = is.get();
//         get_line( is, linebuf );
//         std::vector< char >::const_iterator li = linebuf.begin();
        
        if( c == '>' ) {
            
            
            while( isspace(is.get()) ) {}
            is.unget();
            
            names.push_back(std::string());
            std::string &str = names.back();
        
            while( true ) {
                c = is.get();
                
                if( xisspace(c) || is.eof() ) {
                    break;   
                }
                str.push_back(c); 
            }
            
        
//             std::cout << "name: " << names.back() << std::endl;
            
            data.push_back( std::string() );
            //data_accum = &data.back();
            data_accum = &data.back();
        } else {
      
            if( data_accum == 0 ) {
                throw std::runtime_error( "data_accum == 0. probably bad fasta file\n" );
            }
            
            while( xisspace( is.get() )) {}
         
            is.unget();
            
            while( true ) {
                c = is.get();
                
                if( xisspace(c) || is.eof() ) {
                    break;   
                }
                data_accum->push_back(c); 
                
            }
       }
    }
    
//     std::cout << "size: " << names.size() << " " << data.size() << "\n";
    data.resize( names.size() );
    
}

template<class input, class StateMap>
class inc_fasta {
    input &m_input;
    
    StateMap &m_state_map;
    
public:
    inc_fasta( input &inp, StateMap &state_map ) : m_input(inp), m_state_map(state_map) {
        reset(); 
     
    }
    
    void reset() {
        m_input.clear();
        m_input.seekg(0);   
    }
    
    bool next_seq( std::string &name, std::vector<char> &seq ) {
        while( !m_input.eof() && m_input.get() != '>' ) {}

        if( m_input.eof() ) {
            return false;
        }

        int c;
        while( true ) {
            c = m_input.get();
            
            if( xisspace(c) || m_input.eof() ) {
                break;   
            }
            name.push_back(c); 
        }
        
        
        if( c != '\n' ) {
            while( !m_input.eof() && m_input.get() != '\n' ) {}
        }
            
        while( xisspace( m_input.get() )) {}
        
        m_input.unget();
            
        while( true ) {
            c = m_input.get();
            
            if( m_input.eof() ) {
                break;
            }
            
            if( c == '>' ) {
                m_input.unget();
                break;
            }
                
                
            if( !isspace(c) ) {
                seq.push_back( m_state_map.state_backmap(c)); 
            }
                
        }
        return true;
    }   
    
};

static void write_fasta( std::ostream &os, std::vector<std::string> &names, std::vector<std::string> &data ) {
 
    assert( names.size() == data.size() );
    
    for( size_t i = 0; i < names.size(); i++ ) {
     
        os << ">" << names[i] << "\n";
        os << data[i] << "\n";
    }
        
        
        
}

class scoring_matrix {
public:
    typedef int8_t score_t;
private:
    
    const static size_t MAX_SIZE = 256;
    
    aligned_buffer<score_t> m_cmatrix;
    //std::vector<score_t> m_cmatrix;
    std::vector<score_t> m_matrix;
    std::vector<int> m_backmap;
    std::string m_alphabet;
    size_t addr( int a, int b ) {
        return a + b * MAX_SIZE;
    }
    
    size_t caddr( int a, int b ) {
        return a + b * m_alphabet.size();
    }
    
    void compress() {
        const size_t asize = m_alphabet.size();
        m_cmatrix.resize( asize * asize );
        
        
        for( uint i = 0; i < asize; i++ ) {
            for( uint j = 0; j < asize; j++ ) {
                m_cmatrix.m_ptr[caddr(i,j)] = m_matrix[addr(m_alphabet[i], m_alphabet[j])];
            }
        }
    }
public:
    inline size_t num_states() {
        return m_alphabet.size();
        
    }
    
    inline int get_state( int i ) {
        return m_alphabet[i];
    }

    inline int state_backmap( int s ) {
        return m_backmap[s];
    }

    inline score_t get_score( int a, int b ) {
//         assert( a >= 0 && a < MAX_SIZE && b >= 0 && b < MAX_SIZE );
        return m_matrix[addr(a,b)];
    }
    
    inline score_t *get_slice( int a ) {
        return &m_matrix[a * MAX_SIZE];
    }
    
    inline score_t *get_cslice( int a ) {
        //return &m_cmatrix[a * m_alphabet.size()];
        return m_cmatrix(a * m_alphabet.size());
    }
    
    scoring_matrix( std::istream &is ) 
        : m_matrix( MAX_SIZE * MAX_SIZE ),
            m_backmap(MAX_SIZE)
        
    {
        std::fill(m_backmap.begin(), m_backmap.end(), -1);
        std::vector<char> linebuf;
        
        bool have_firstline = false;
        size_t mline = 0;
        while( is.good() && !is.eof() ) {
            
            char c;
            
           
            do {
                c = is.get();
            } while( std::isspace(c) );
        
            if( is.eof() ) {
                break;
            }
            
//            get_line( is, linebuf );
                
            if( c == '#' ) {
                while( !is.eof() && is.get() != '\n' ) {}
                continue;
            }
            
  
            
            if( !have_firstline ) {
                while( c != '\n' ) {
                    m_backmap[c] = m_alphabet.size();
                    m_alphabet.push_back(c);
                    
                    do {
                        c = is.get();
                        
                        if( c == '\n' ) {
                            break;
                        }
                    } while( std::isspace(c) );
                }
            
//                 std::cout << "alphabet: " << m_alphabet << "\n";
                have_firstline = true;
            } else {
                is.unget();
//                 char c = is.get();
//                 std::cout << "c: " << c << "\n";
                
                std::string lname;
                is >> lname;
                
//                 std::cout << "name: " << lname << "\n";
                
                char lnc = lname[0];
                
                for( size_t i = 0; i < m_alphabet.size(); i++ ) {
                 
                    int score;
                    
                    is >> score;
                    
//                     std::cout << "score: " << int(lnc) << " " << int(m_alphabet[i]) << " " << score << "\n";
                    m_matrix[addr( lnc, m_alphabet[i] )] = score;
                    m_matrix[addr( m_alphabet[i], lnc )] = score;
                }
            }
        }
     
        compress();
//         for( size_t i = 0; i < m_alphabet.size(); i++ ) {
//             for( size_t j = 0; j < m_alphabet.size(); j++ ) {
//                 std::cout << int(get_score( m_alphabet[i], m_alphabet[j] )) << " ";
//             }
//             std::cout << "\n";
//             
//             
//         }
     
    }
};
