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


#include "PsiData.hh"
#include "Comm.hh"
#include "Mat.hh"
#include <string.h>
#include <stdio.h>

/*
    writePsiToFile

    Writes psi to a file in parallel.

    Data Format:
    char[32]: "Tycho 2 Psi Output" (zeros for any trailing characters
    uint64_t: version of file format
    uint64_t: number of cells
    uint64_t: number of angles
    uint64_t: number of energy groups
    CellData[]: array of cell data

    CellData Format:
    double[]: psi(:, :, :, global cell index)
*/
void writePsiToFile(const std::string &filename,
                    const PsiData &psi)
{
  FILE *ofp;

  ofp = fopen(filename.c_str(),"w");

  for (size_t cell = 0; cell < g_nCells; cell++) {
    for (UINT a = 0; a < g_nAngles; a++) {
      for (UINT v = 0; v < g_nVrtxPerCell; v++) {
        for (UINT g = 0; g < g_nGroups; g++) {
            fprintf(ofp,"psi(%lu, %lu, %lu, %d): %4f \n", psi(g,v,a,cell));    
        }}}}
  fclose(ofp); 

#if 0
    Comm::File file;
    char outputName[32] = {
        'T', 'y', 'c', 'h', 'o', ' ', '2', ' ', 'P', 's', 'i', ' ', 
        'O', 'u', 't', 'p', 'u', 't', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    uint64_t restOfHeader[4];


    // Fill in rest of header
    // version, number of cells, number of angles, and number of group
    restOfHeader[0] = 1;
    restOfHeader[1] = g_nCells;
    Comm::gsum(restOfHeader[1]);
    restOfHeader[2] = g_nAngles;
    restOfHeader[3] = g_nGroups;


    // Open file
    std::string psiFile(filename+".psi");
    Comm::openFileForWrite(psiFile, file);

    // If rank 0, write header data
    if (Comm::rank() == 0) {
        double header[8];
        memcpy(header, outputName, 4 * sizeof(double));
        memcpy(&header[4], restOfHeader, 4 * sizeof(double));
        Comm::writeDoublesAt(file, 0, header, 8);
    }


    // Write data one cell at a time
    Mat3<double> cellData(g_nGroups, g_nVrtxPerCell, g_nAngles);
    for (size_t cell = 0; cell < g_nCells; cell++) {
        int dataSize = g_nAngles * g_nGroups * g_nVrtxPerCell;
        uint64_t globalCell = g_tychoMesh->getLGCell(cell);
        uint64_t offset = 8 + globalCell * dataSize;
        
        for (UINT a = 0; a < g_nAngles; a++) {
        for (UINT v = 0; v < g_nVrtxPerCell; v++) {
        for (UINT g = 0; g < g_nGroups; g++) {
            cellData(g,v,a) = psi(g,v,a,cell);    
        }}}
        double *data = &cellData[0];
        Comm::writeDoublesAt(file, offset, data, dataSize);
    }


    // Close file
    Comm::closeFile(file);
#endif
}




void writePhiToFile(const std::string &filename,
                    const PhiData &phi)
{
  FILE *ofp;
  ofp = fopen(filename.c_str(),"w");

  for (size_t cell = 0; cell < g_nCells; cell++) {
      for (UINT v = 0; v < g_nVrtxPerCell; v++) {
        for (UINT g = 0; g < g_nGroups; g++) {
            fprintf(ofp,"phi(%lu, %lu, %d): %4f \n", phi(g,v,cell));    
        }}}
  fclose(ofp); 


#if 0
    Comm::File file;
    char outputName[32] = {
        'T', 'y', 'c', 'h', 'o', ' ', '2', ' ', 'P', 'h', 'i', ' ', 
        'O', 'u', 't', 'p', 'u', 't', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    uint64_t restOfHeader[4];


    // Fill in rest of header
    // version, number of cells, number of angles, and number of group
    restOfHeader[0] = 1;
    restOfHeader[1] = g_nCells;
    Comm::gsum(restOfHeader[1]);
    restOfHeader[2] = g_nAngles;
    restOfHeader[3] = g_nGroups;


    // Open file
    std::string phiFile(filename+".phi");
    Comm::openFileForWrite(filename, file);
    

    // If rank 0, write header data
    if (Comm::rank() == 0) {
        double header[8];
        memcpy(header, outputName, 4 * sizeof(double));
        memcpy(&header[4], restOfHeader, 4 * sizeof(double));
        Comm::writeDoublesAt(file, 0, header, 8);
    }


    // Write data one cell at a time
    Mat2<double> cellData(g_nGroups, g_nVrtxPerCell);
    for (size_t cell = 0; cell < g_nCells; cell++) {
        int dataSize = g_nGroups * g_nVrtxPerCell;
        uint64_t globalCell = g_tychoMesh->getLGCell(cell);
        uint64_t offset = 8 + globalCell * dataSize;
        
        for (UINT v = 0; v < g_nVrtxPerCell; v++) {
        for (UINT g = 0; g < g_nGroups; g++) {
            cellData(g,v) = phi(g,v,cell);    
        }}
        double *data = &cellData[0];
        Comm::writeDoublesAt(file, offset, data, dataSize);
    }


    // Close file
    Comm::closeFile(file);
#endif
}
