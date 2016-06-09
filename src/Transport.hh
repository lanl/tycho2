//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   transport/Transport.hh
 * \author Shawn Pautz
 * \date   Fri Jan 14 16:30:59 2000
 * \brief  Transport class header file.
 */
//---------------------------------------------------------------------------//
// $Id: Transport.hh,v 1.12 2000/10/17 16:36:45 pautz Exp $
//---------------------------------------------------------------------------//

#ifndef __transport_Transport_hh__
#define __transport_Transport_hh__

#include "Mat.hh"
#include "PsiData.hh"
#include "Typedef.hh"
#include <vector>


namespace Transport 
{
    void solve(const UINT cell, const UINT angle, 
               const double sigmaTotal,
               const Mat3<double> &localPsiBound, 
               const Mat2<double> &localSource,
               Mat2<double> &localPsi);
} // End namespace Transport

#endif

