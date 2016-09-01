/*
    PsiData.hh
    
    Class for data for psi.  Want the ability to change order of data
    data structures.
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


#ifndef __PSI_DATA_HH__
#define __PSI_DATA_HH__

#include "Assert.hh"
#include "Typedef.hh"
#include "Global.hh"
#include "Quadrature.hh"
#include "TychoMesh.hh"
#include <stddef.h>


/*
    PsiData

    g = group
    v = vertex
    a = angle
    c = cell
*/
class PsiData {
public:
    
    // Accessors
    double& operator()(size_t g, size_t v, size_t a, size_t c) 
    {
        return c_data[index(g,v,a,c)];
    }
    
    const double& operator()(size_t g, size_t v, size_t a, size_t c) const 
    {
        return c_data[index(g,v,a,c)];
    }
    
    double& operator[](size_t i)
    {
        check(i);
        return c_data[i];
    }

    const double& operator[](size_t i) const
    {
        check(i);
        return c_data[i];
    }


    // Size of data structure
    size_t size() const { return c_ng * c_nv * c_na * c_nc; }


    // Constructor
    PsiData(double initValue = 0.0)
    {
        c_ng = g_nGroups;
        c_nv = g_nVrtxPerCell;
        c_na = g_quadrature->getNumAngles();
        c_nc = g_spTychoMesh->getNCells();
        c_data = new double[size()];
        setToValue(initValue);
    }
    
    // Copy constructor
    PsiData(const PsiData &other)
    {
        c_ng = g_nGroups;
        c_nv = g_nVrtxPerCell;
        c_na = g_quadrature->getNumAngles();
        c_nc = g_spTychoMesh->getNCells();
        c_data = new double[size()];
        for (size_t i = 0; i < size(); i++) {
            c_data[i] = other.c_data[i];
        }
    }

    // Assignment Operator
    PsiData & operator= (const PsiData &other)
    {
        if (this != &other) {
            if (c_data != NULL) {
               delete[] c_data;
            }

            c_ng = g_nGroups;
            c_nv = g_nVrtxPerCell;
            c_na = g_quadrature->getNumAngles();
            c_nc = g_spTychoMesh->getNCells();
            c_data = new double[size()];
            for (size_t i = 0; i < size(); i++) {
                c_data[i] = other.c_data[i];
            }
        }    
    }


    //Destructor
    ~PsiData()
    {
        if (c_data != NULL) {
            delete[] c_data;
            c_data = NULL;
        }
    }
    
    
    // Set constant value
    void setToValue(double value)
    {
        if(c_data == NULL)
            return;
        
        for(size_t i = 0; i < size(); i++) {
            c_data[i] = value;
        }
    }


// Private    
private:
    size_t c_ng, c_nv, c_na, c_nc;
    double *c_data;


    // Compute the offset into the data array.
    size_t index(size_t g, size_t v, size_t a, size_t c) const
    {
        Assert(g < c_ng);
        Assert(v < c_nv);
        Assert(a < c_na);
        Assert(c < c_nc);
        
        return ((c * c_na + a) * c_nv + v) * c_ng + g;
    }


    // Make sure a bare integer index is within the appropriate range.
    void check(size_t i) const
    {
        UNUSED_VARIABLE(i);     // Needed when Assert does nothing
        Assert( i < size() );
    }
};


/*
    PsiBoundData

    g = group
    v = vertex
    a = angle
    s = side
*/
class PsiBoundData {
public:
    
    // Accessors
    double& operator()(size_t g, size_t v, size_t a, size_t s) 
    {
        return c_data[index(g,v,a,s)];
    }
    
    const double& operator()(size_t g, size_t v, size_t a, size_t s) const 
    {
        return c_data[index(g,v,a,s)];
    }
    
    double& operator[](size_t i)
    {
        check(i);
        return c_data[i];
    }

    const double& operator[](size_t i) const
    {
        check(i);
        return c_data[i];
    }


    // Size of data structure
    size_t size() const { return c_ng * c_nv * c_na * c_ns; }


    // Constructor
    PsiBoundData(double initValue = 0.0)
    {
        c_ng = g_nGroups;
        c_nv = g_nVrtxPerFace;
        c_na = g_quadrature->getNumAngles();
        c_ns = g_spTychoMesh->getNSides();
        c_data = new double[size()];
        setToValue(initValue);
    }
    
    // Copy constructor
    PsiBoundData(const PsiBoundData &other)
    {
        c_ng = g_nGroups;
        c_nv = g_nVrtxPerFace;
        c_na = g_quadrature->getNumAngles();
        c_ns = g_spTychoMesh->getNSides();
        c_data = new double[size()];
        for (size_t i = 0; i < size(); i++) {
            c_data[i] = other.c_data[i];
        }
    }

    // Assignment Operator
    PsiBoundData & operator= (const PsiBoundData &other)
    {
        if (this != &other) {
            if (c_data != NULL) {
               delete[] c_data;
            }

            c_ng = g_nGroups;
            c_nv = g_nVrtxPerFace;
            c_na = g_quadrature->getNumAngles();
            c_ns = g_spTychoMesh->getNSides();
            c_data = new double[size()];
            for (size_t i = 0; i < size(); i++) {
                c_data[i] = other.c_data[i];
            }
        }    
    }


    //Destructor
    ~PsiBoundData()
    {
        if (c_data != NULL) {
            delete[] c_data;
            c_data = NULL;
        }
    }
    
    
    // Set constant value
    void setToValue(double value)
    {
        if(c_data == NULL)
            return;
        
        for(size_t i = 0; i < size(); i++) {
            c_data[i] = value;
        }
    }


// Private    
private:
    size_t c_ng, c_nv, c_na, c_ns;
    double *c_data;


    // Compute the offset into the data array.
    size_t index(size_t g, size_t v, size_t a, size_t s) const
    {
        Assert(g < c_ng);
        Assert(v < c_nv);
        Assert(a < c_na);
        Assert(s < c_ns);
        
        return ((s * c_na + a) * c_nv + v) * c_ng + g;
    }


    // Make sure a bare integer index is within the appropriate range.
    void check(size_t i) const
    {
        UNUSED_VARIABLE(i);     // Needed when Assert does nothing
        Assert( i < size() );
    }
};







class PhiData {
public:
    
    // Accessors
    double& operator()( size_t i, size_t j, size_t k ) 
        { return c_data[index(i,j,k)]; }
    
    const double& operator()( size_t i, size_t j, size_t k ) const 
        { return c_data[index(i,j,k)]; }
    
    double& operator[]( size_t i ) { check(i); return c_data[i]; }
    const double& operator[]( size_t i ) const { check(i); return c_data[i]; }


    // Size of Mat
    size_t size() const { return c_len1 * c_len2 * c_len3; }


    // Constructor
    PhiData( size_t len1, size_t len2, size_t len3, double initValue = 0.0 )
    {
        c_len1 = len1;
        c_len2 = len2;
        c_len3 = len3;
        c_data = new double[size()];
        setToValue(initValue);
    }

    //Copy constructor
    PhiData(const PhiData &other)
    {
        c_len1=other.c_len1;
        c_len2=other.c_len2;
        c_len3=other.c_len3;
        c_data = new double[size()];
        for (size_t i=0; i<size(); i++)
            {c_data[i] = other.c_data[i];}
    }

    //Assignment Operator
    PhiData & operator= (const PhiData &other) //= delete;
    {
        if (this != &other)
        {
            if (c_data != NULL) {
               delete[] c_data;
            }
            c_len1=other.c_len1;
            c_len2=other.c_len2;
            c_len3=other.c_len3;
            c_data = new double[size()];
            for (size_t i=0; i<size(); i++)
                {c_data[i] = other.c_data[i];}

        }    
    }

    // Destructor
    ~PhiData()
    {
        if (c_data != NULL) {
            delete[] c_data;
            c_data = NULL;
        }
    }
    
    
    void setToValue(double value)
    {
        if(c_data == NULL)
            return;
        
        for(size_t i = 0; i < size(); i++) {
            c_data[i] = value;
        }
    }


// Private    
private:
    size_t c_len1, c_len2, c_len3;
    double *c_data;

    // Compute the offset into the data array.
    size_t index( size_t i, size_t j, size_t k ) const
    {
        //Assert( i >= 0);
        Assert( i < c_len1 );
        //Assert( j >= 0);
        Assert( j < c_len2 );
        //Assert( k >= 0);
        Assert( k < c_len3 );
        
        return c_len1 * (c_len2 * k + j) + i;
    }

    // Make sure a bare integer index is within the appropriate range.
    void check( size_t i ) const
    {
        UNUSED_VARIABLE(i);
        //Assert( i >= 0 );
        Assert( i < size() );
    }
};


#endif
