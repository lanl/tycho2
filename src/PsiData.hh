//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PsiData.hh
 * \author Kris Garrett
 * \date   December 2015
 * \brief  Class for data for psi.  Want the ability to change order of data
 *         data structures.
 */
//---------------------------------------------------------------------------//
// $Id: Mat.hh,v 1.17 2003/03/14 19:55:03 tme Exp $
//---------------------------------------------------------------------------//

#ifndef __PSI_DATA_HH__
#define __PSI_DATA_HH__

#include "Assert.hh"
#include "Typedef.hh"
#include <stddef.h>

class PsiData {
public:
    
    // Accessors
    double& operator()( size_t i, size_t j, size_t k, size_t l ) 
        { return c_data[index(i,j,k,l)]; }
    
    const double& operator()( size_t i, size_t j, size_t k, size_t l ) const 
        { return c_data[index(i,j,k,l)]; }
    
    double& operator[]( size_t i ) { check(i); return c_data[i]; }
    const double& operator[]( size_t i ) const { check(i); return c_data[i]; }


    // Size of Mat
    size_t size() const { return c_len1 * c_len2 * c_len3 * c_len4; }


    // Constructor
    PsiData( size_t len1, size_t len2, size_t len3, size_t len4, double initValue = 0.0 )
    {
        c_len1 = len1;
        c_len2 = len2;
        c_len3 = len3;
        c_len4 = len4;
        c_data = new double[size()];
        setToValue(initValue);
    }


    // Destructor
    ~PsiData()
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
    size_t c_len1, c_len2, c_len3, c_len4;
    double *c_data;

    // Compute the offset into the data array.
    size_t index( size_t i, size_t j, size_t k, size_t l ) const
    {
        //Assert( i >= 0);
        Assert( i < c_len1 );
        //Assert( j >= 0);
        Assert( j < c_len2 );
        //Assert( k >= 0);
        Assert( k < c_len3 );
        //Assert( l >= 0);
        Assert( l < c_len4 );
        
        return c_len1 * (c_len2 * (c_len3 * l + k) + j) + i;
    }

    // Make sure a bare integer index is within the appropriate range.
    void check( size_t i ) const
    {
        UNUSED_VARIABLE(i);
        //Assert( i >= 0 );
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
