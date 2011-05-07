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
#ifndef __aligned_buffer_h
#define __aligned_buffer_h


#include <cstdlib>
#include <cstddef>
#include <stdexcept>

template<class T>
struct aligned_buffer {
    T* m_ptr;
    size_t m_size;
    const static size_t align = 32;
    
    aligned_buffer() : m_ptr(0), m_size(0) {}
    
    aligned_buffer( size_t size ) : m_ptr(0), m_size(0) {
         
        resize( size );
                
    }
    
    void resize( size_t size ) {
    
        free( m_ptr );

        m_size = size;
        int ret = posix_memalign( (void**)&m_ptr, align, byte_size() );
        
        if( ret != 0 ) {
            throw std::runtime_error( "posix_memalign failed" );
        }

        
    }
    
    size_t size() {
        return m_size;
    }
    
    size_t byte_size() {
        return m_size * sizeof(T);
    }
    
    ~aligned_buffer() {
        free( m_ptr );
    }
    
    T *begin() {
        return m_ptr;
    }
    
    T* end() {
        return begin() + m_size;
    }
    
    inline T* operator() (ptrdiff_t o) {
        return begin() + o;
    }
    
private:
    aligned_buffer( const aligned_buffer &other ) {}
    aligned_buffer &operator=( const aligned_buffer &other ) {}
    
};

#endif