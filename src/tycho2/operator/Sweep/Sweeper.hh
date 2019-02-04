/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.
*/

#include "PsiData.hh"
#include "SweeperAbstract.hh"

class Sweeper : public SweeperAbstract
{
public:
    Sweeper();
    void sweep(PsiData &psi, const PsiData &source, bool zeroPsiBound);
    void solve();

private:

};
