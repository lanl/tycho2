/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.
*/

#ifndef __SWEEPER_ABSTRACT_HH__
#define __SWEEPER_ABSTRACT_HH__

#include "PsiData.hh"
#include <string>


class SweeperAbstract
{
public:
    virtual
    void sweep(PsiData &psi, const PsiData &source, 
               bool zeroPsiBound = false) = 0;

    virtual
    void solve() = 0;

    void writePsiToFile(std::string &filename)
    {
        c_psi.writeToFile(filename);
    }

    PsiData& getPsi()
    {
        return c_psi;
    }

protected:
    PsiData c_psi;
    PsiData c_source;
};

#endif
