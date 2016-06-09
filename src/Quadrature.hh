/*
    Quadrature.hh
*/

#ifndef __QUADRATURE_HH__
#define __QUADRATURE_HH__

#include "Typedef.hh"
#include <vector>


class Quadrature
{
public:
    Quadrature(const UINT snOrder);
    
    double getMu(const UINT angle) const;
    double getEta(const UINT angle) const;
    double getXi(const UINT angle) const;
    double getWt(const UINT angle) const;
    std::vector<double> getOmega(const UINT angle) const;
    UINT getNumAngles() const;
    
private:
    UINT c_numAngles;
    std::vector<double> c_xi;
    std::vector<double> c_eta;
    std::vector<double> c_mu;
    std::vector<double> c_w;
};


#endif

