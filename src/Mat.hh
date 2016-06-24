/*
    Mat.hh
    
    Implements 1, 2, and 3 dimensional arrays.
*/

/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced 
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for 
the U.S. Department of Energy. The U.S. Government has rights to use, 
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified 
to produce derivative works, such modified software should be clearly marked, 
so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
1.      Redistributions of source code must retain the above copyright notice, 
        this list of conditions and the following disclaimer.
2.      Redistributions in binary form must reproduce the above copyright 
        notice, this list of conditions and the following disclaimer in the 
        documentation and/or other materials provided with the distribution.
3.      Neither the name of Los Alamos National Security, LLC, Los Alamos 
        National Laboratory, LANL, the U.S. Government, nor the names of its 
        contributors may be used to endorse or promote products derived from 
        this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND 
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL 
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef __MAT_HH__
#define __MAT_HH__

#include "Assert.hh"
#include "Typedef.hh"
#include <algorithm>
#include <memory>
#include <map>
#include <vector>


/*
    Mat1
    
    Implements 1D Array
*/
template< class T >
class Mat1 {
private:
    
    // Private data
    size_t c_xlen;
    T *c_v;
    
    
    // Private functions
    size_t index( size_t i ) const
    {
        Assert( i >= 0 );
        Assert( i < c_xlen );
        return i;
    }
    
    void detach()
    {
        if (c_v != NULL) {
            delete[] c_v;
            c_v = NULL;
        }
    }
    

public:
    
    // Accessors
    T&       operator[]( size_t i )       { return c_v[ index(i) ]; }
    const T& operator[]( size_t i ) const { return c_v[ index(i) ]; }

    T&       operator()( size_t i )       { return c_v[ index(i) ]; }
    const T& operator()( size_t i ) const { return c_v[ index(i) ]; }


    // Size of Mat
    size_t size() const { return c_xlen; }
    
    
    // Constructors
    Mat1()
    {
        c_xlen = 0;
        c_v = NULL;
    }
    
    Mat1( size_t xmax, const T& t = T() )
    {
        c_xlen = xmax;
        c_v = new T[size()];
        std::uninitialized_fill( c_v, c_v + size(), t );
    }

    Mat1( const Mat1<T>& m )
    {
        c_xlen = m.c_xlen;
        c_v = new T[size()];
        std::uninitialized_copy( m.c_v, m.c_v + m.size(), c_v );
    }
    
    
    // Destructor
    ~Mat1()
    {
        detach();
    }
    
    
    // Assignment operators
    Mat1& operator=( const T& t )
    {
        Assert(c_v != NULL);
        std::fill( c_v, c_v + size(), t );
        return *this;
    }
    
    Mat1& operator=( const Mat1& m )
    {
        if (this == &m)
            return *this;

        if ( m.c_xlen != c_xlen ) {
            detach();
            c_xlen = m.c_xlen;
            c_v = new T[size()];
            std::uninitialized_copy( m.c_v, m.c_v + m.size(), c_v );
        }
        else {
            if (c_v != NULL)
                std::copy( m.c_v, m.c_v + m.size(), c_v );
        }

        return *this;
    }

    
    // Arithmetic
    Mat1& operator+=( const Mat1<T>& m )
    {
        Assert(c_xlen == m.c_xlen);
        for(size_t i = 0; i < size(); i++) {
            c_v[i] += m.c_v[i];
        }
        return *this;
    }
    

    // Utility support
    void redim( size_t nxmax, const T& t = T() )
    {
        if (c_v != NULL) {
            detach();
        }
        c_xlen = nxmax;
        c_v = new T[size()];
        std::uninitialized_fill( c_v, c_v + size(), t );
    }
};


/*
    Mat2
    
    Implements 2D Array
*/

template< class T >
class Mat2 {
private:
    size_t c_xlen, c_ylen;
    T *c_v;

    // Compute the offset into the data array, of the i,j th element.
    size_t index( size_t i, size_t j ) const
    {
        Assert( i >= 0 );
        Assert( i < c_xlen );
        Assert( j >= 0 );
        Assert( j < c_ylen );

        return c_xlen * j + i;
    }

    // Make sure a bare integer index is within the appropriate range.
    void check( size_t i ) const
    {
        UNUSED_VARIABLE(i);
        Assert( i >= 0 );
        Assert( i < size() );
    }

    void detach()
    {
        if (c_v != NULL) {
            delete[] c_v;
            c_v = NULL;
        }
    }

public:

    // Accessors
    T&       operator()( size_t i, size_t j )       { return c_v[ index(i,j) ]; }
    const T& operator()( size_t i, size_t j ) const { return c_v[ index(i,j) ]; }

    T& operator[]( size_t i ) { check(i); return c_v[i]; }
    const T& operator[]( size_t i ) const { check(i); return c_v[i]; }
    
    
    // Size of Mat
    size_t size() const { return c_xlen * c_ylen; }
    
    
    // Constructors
    Mat2()
    {
        c_xlen = 0;
        c_ylen = 0;
        c_v = NULL;
    }

    Mat2( size_t xmax, size_t ymax, const T& t = T() )
    {
        c_xlen = xmax;
        c_ylen = ymax;
        c_v = new T[size()];
        std::uninitialized_fill( c_v, c_v + size(), t );
    }
    
    Mat2( size_t xmax, size_t ymax, const T *data )
    {
        c_xlen = xmax;
        c_ylen = ymax;
        c_v = new T[size()];
        for (size_t i = 0; i < xmax * ymax; i++) 
            c_v[i] = data[i];
    }

    Mat2( const Mat2<T>& m )
    {
        c_xlen = m.c_xlen;
        c_ylen = m.c_ylen;
        c_v = new T[size()];
        std::uninitialized_copy( m.c_v, m.c_v + m.size(), c_v );
    }


    // Destructor
    ~Mat2()
    {
        detach();
    }


    // Assignment operators
    Mat2& operator=( const T& t )
    {
        std::fill( c_v, c_v + size(), t );
        return *this;
    }

    Mat2& operator=( const Mat2& m )
    {
        if (this == &m) return *this;

        if ( m.c_xlen != c_xlen ||
             m.c_ylen != c_ylen ) {
            detach();
            c_xlen = m.c_xlen;
            c_ylen = m.c_ylen;
            c_v = new T[size()];
            std::uninitialized_copy( m.c_v, m.c_v + m.size(), c_v );
        }
        else {
            if (c_v)
            std::copy( m.c_v, m.c_v + m.size(), c_v );
        }

        return *this;
    }

    
    // Arithmetic
    Mat2& operator+=( const Mat2<T>& m )
    {
        Assert( c_xlen == m.c_xlen );
        Assert( c_ylen == m.c_ylen );
        
        for(size_t i = 0; i < size(); i++) {
            c_v[i] += m.c_v[i];
        }
        return *this;
    }
    
    
    // Utility Support
    void redim( size_t nxmax, size_t nymax, const T& t = T() )
    {
        if (c_v != NULL) {
            detach();
        }
        c_xlen = nxmax;
        c_ylen = nymax;
        c_v = new T[size()];
        std::uninitialized_fill( c_v, c_v + size(), t );
    }
};


/*
    Mat3
    
    Implements 3D Array
*/
template< class T >
class Mat3 {
private:
    size_t c_xlen, c_ylen, c_zlen;
    T *c_v;

    // Compute the offset into the data array, of the i,j th element.
    size_t index( size_t i, size_t j, size_t k ) const
    {
        Assert( i >= 0 );
        Assert( i < c_xlen );
        Assert( j >= 0 );
        Assert( j < c_ylen );
        Assert( k >= 0 );
        Assert( k < c_zlen );

        return c_xlen * (k * c_ylen + j) + i;
    }

    // Make sure a bare integer index is within the appropriate range.
    void check( size_t i ) const
    {
        Assert( i >= 0 );
        Assert( i < size() );
    }

    void detach()
    {
        if (c_v != NULL) {
            delete[] c_v;
            c_v = NULL;
        }
    }
    

public:

    // Accessors
    T&       operator()( size_t i, size_t j, size_t k )       { return c_v[ index(i,j,k) ]; }
    const T& operator()( size_t i, size_t j, size_t k ) const { return c_v[ index(i,j,k) ]; }

    T& operator[]( size_t i ) { check(i); return c_v[i]; }
    const T& operator[]( size_t i ) const { check(i); return c_v[i]; }


    // Size of Mat
    size_t size() const { return c_xlen * c_ylen * c_zlen; }


    // Constructors
    Mat3()
    {
        c_xlen = 0;
        c_ylen = 0;
        c_zlen = 0;
        c_v = NULL;
    }

    Mat3( size_t xmax, size_t ymax, size_t zmax, const T& t = T() )
    {
        c_xlen = xmax;
        c_ylen = ymax;
        c_zlen = zmax;
        c_v = new T[size()];
        std::uninitialized_fill( c_v, c_v + size(), t );
    }

    Mat3( const Mat3<T>& m )
    {
        c_xlen = m.c_xlen;
        c_ylen = m.c_ylen;
        c_zlen = m.c_zlen;
        c_v = new T[size()];
        std::uninitialized_copy( m.c_v, m.c_v + m.size(), c_v );
    }


    // Destructor
    ~Mat3()
    {
        detach();
    }


    // Assignment operators
    Mat3& operator=( const T& t )
    {
        std::fill( c_v, c_v + size(), t );
        return *this;
    }

    Mat3& operator=( const Mat3& m )
    {
        if (this == &m) return *this;

        if ( m.c_xlen != c_xlen ||
             m.c_ylen != c_ylen ||
             m.c_zlen != c_zlen ) {
            detach();
            c_xlen = m.c_xlen;
            c_ylen = m.c_ylen;
            c_zlen = m.c_zlen;
            c_v = new T[size()];
            std::uninitialized_copy( m.c_v, m.c_v + m.size(), c_v );
        }
        else {
            if (c_v)
            std::copy( m.c_v, m.c_v + m.size(), c_v );
        }

        return *this;
    }

    
    // Arithmetic
    Mat3& operator+=( const Mat3<T>& m )
    {
        Assert( c_xlen == m.c_xlen );
        Assert( c_ylen == m.c_ylen );
        Assert( c_zlen == m.c_zlen );
        
        for(size_t i = 0; i < size(); i++) {
            c_v[i] += m.c_v[i];
        }
        return *this;
    }


    // Utility support
    void redim( size_t nxmax, size_t nymax, size_t nzmax, const T& t = T() )
    {
        if (c_v != NULL) {
            detach();
        }
        c_xlen = nxmax;
        c_ylen = nymax;
        c_zlen = nzmax;
        c_v = new T[size()];
        std::uninitialized_fill( c_v, c_v + size(), t );
    }
};



// extern template class Mat1<double>;
// extern template class Mat1<bool>;
// extern template class Mat1<UINT>;
// 
// extern template class Mat2<double>;
// extern template class Mat2<bool>;
// extern template class Mat2<UINT>;
// 
// extern template class Mat3<double>;
// extern template class Mat3<bool>;
// extern template class Mat3<UINT>;
// 
// extern template class std::map<UINT,UINT>;
// extern template class std::vector<double>;
// extern template class std::vector<int>;
// extern template class std::vector<UINT>;


#endif

