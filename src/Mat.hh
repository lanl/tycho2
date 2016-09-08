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
#include <stddef.h>


/*
    Mat1
    
    Implements 1D Array
*/
template< class T >
class Mat1
{
private:
    
    // Data members
    size_t c_xlen;
    T *c_v;
    
    
    // Check size before returning
    size_t index(size_t i) const
    {
        Assert(i < c_xlen);
        return i;
    }
    

    // Delete data
    void detach()
    {
        if (c_v != NULL) {
            delete[] c_v;
            c_v = NULL;
        }
    }
    

public:
    
    // Accessors
    T& operator[](size_t i)
    {
        return c_v[index(i)];
    }
    const T& operator[](size_t i) const
    {
        return c_v[index(i)];
    }

    T& operator()(size_t i)
    {
        return c_v[index(i)];
    }
    const T& operator()(size_t i) const
    {
        return c_v[index(i)];
    }


    // Size of Mat
    size_t size() const
    {
        return c_xlen;
    }
    
    
    // Constructors
    Mat1()
    {
        c_xlen = 0;
        c_v = NULL;
    }
    
    Mat1(size_t xmax)
    {
        c_xlen = xmax;
        c_v = new T[size()];

        for (size_t i = 0; i < size(); i++) {
            c_v[i] = T();
        }
    }

    
    // Destructor
    ~Mat1()
    {
        detach();
    }
    
    
    // Don't allow copy or assignment operators
    Mat1(const Mat1<T> &m) = delete;
    Mat1& operator=(const Mat1 &m) = delete;
    
    
    // Resize matrix
    void resize(size_t nxmax)
    {
        detach();
        c_xlen = nxmax;
        c_v = new T[size()];
        
        for (size_t i = 0; i < size(); i++) {
            c_v[i] = T();
        }
    }
};


/*
    Mat2
    
    Implements 2D Array
*/

template< class T >
class Mat2
{
private:
    
    // Data members
    size_t c_xlen, c_ylen;
    T *c_v;


    // Compute the offset into the data array, of the (i,j) element.
    size_t index(size_t i, size_t j) const
    {
        Assert(i < c_xlen);
        Assert(j < c_ylen);

        return j * c_xlen + i;
    }


    // Delete data
    void detach()
    {
        if (c_v != NULL) {
            delete[] c_v;
            c_v = NULL;
        }
    }

public:

    // Accessors
    T& operator()(size_t i, size_t j)
    {
        return c_v[index(i, j)];
    }
    const T& operator()(size_t i, size_t j) const
    {
        return c_v[index(i, j)];
    }

    T& operator[](size_t i)
    {
        Assert(i < size());
        return c_v[i];
    }
    const T& operator[](size_t i) const
    {
        Assert(i < size());
        return c_v[i];
    }
    
    
    // Size of Mat
    size_t size() const
    {
        return c_xlen * c_ylen;
    }
    
    
    // Constructors
    Mat2()
    {
        c_xlen = 0;
        c_ylen = 0;
        c_v = NULL;
    }

    Mat2(size_t xmax, size_t ymax)
    {
        c_xlen = xmax;
        c_ylen = ymax;
        c_v = new T[size()];
        
        for (size_t i = 0; i < size(); i++) {
            c_v[i] = T();
        }
    }
    
    
    // Destructor
    ~Mat2()
    {
        detach();
    }


    // Don't allow copy or assignment operators
    Mat2(const Mat2<T> &m) = delete;
    Mat2& operator=(const Mat2 &m) = delete;
    
    
    // Set all values to a constant
    void setAll(const T &t)
    {
        for (size_t i = 0; i < size(); i++) {
            c_v[i] = t;
        }
    }

    
    // Set raw data pointer
    void setData(const T *data)
    {
        for (size_t i = 0; i < size(); i++) {
            c_v[i] = data[i];
        }
    }
    

    // Resize matrix
    void resize( size_t nxmax, size_t nymax)
    {
        detach();
        c_xlen = nxmax;
        c_ylen = nymax;
        c_v = new T[size()];

        for (size_t i = 0; i < size(); i++) {
            c_v[i] = T();
        }
    }
};


/*
    Mat3
    
    Implements 3D Array
*/
template< class T >
class Mat3
{
private:
    
    // Data members
    size_t c_xlen, c_ylen, c_zlen;
    T *c_v;


    // Compute the offset into the data array, of the (i,j,k) element.
    size_t index(size_t i, size_t j, size_t k) const
    {
        Assert(i < c_xlen);
        Assert(j < c_ylen);
        Assert(k < c_zlen);

        return (k * c_ylen + j) * c_xlen + i;
    }

    
    // Delete data
    void detach()
    {
        if (c_v != NULL) {
            delete[] c_v;
            c_v = NULL;
        }
    }
    

public:

    // Accessors
    T& operator()(size_t i, size_t j, size_t k)
    {
        return c_v[index(i, j, k)];
    }
    const T& operator()(size_t i, size_t j, size_t k) const
    {
        return c_v[index(i, j, k)];
    }

    T& operator[](size_t i)
    {
        Assert(i < size());
        return c_v[i];
    }
    const T& operator[](size_t i) const
    {
        Assert(i < size());
        return c_v[i];
    }


    // Size of Mat
    size_t size() const
    {
        return c_xlen * c_ylen * c_zlen;
    }


    // Constructors
    Mat3()
    {
        c_xlen = 0;
        c_ylen = 0;
        c_zlen = 0;
        c_v = NULL;
    }

    Mat3(size_t xmax, size_t ymax, size_t zmax)
    {
        c_xlen = xmax;
        c_ylen = ymax;
        c_zlen = zmax;
        c_v = new T[size()];

        for (size_t i = 0; i < size(); i++) {
            c_v[i] = T();
        }
    }

    
    // Destructor
    ~Mat3()
    {
        detach();
    }


    // Don't allow copy or assignment operators
    Mat3(const Mat3<T> &m) = delete;
    Mat3& operator=(const Mat3 &m) = delete;
    

    // Resize matrix
    void resize(size_t nxmax, size_t nymax, size_t nzmax)
    {
        detach();
        c_xlen = nxmax;
        c_ylen = nymax;
        c_zlen = nzmax;
        c_v = new T[size()];
        
        for (size_t i = 0; i < size(); i++) {
            c_v[i] = T();
        }
    }
};


#endif

