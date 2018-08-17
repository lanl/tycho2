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
#include <string.h>


/*
    writeToFile

    Writes psi to a file in parallel.

    Data Format:
    char[32]: "Tycho 2 Psi Output" (zeros for any trailing characters
    uint64_t: version of file format
    uint64_t: number of cells
    uint64_t: number of angles
    uint64_t: number of energy groups
    CellData[]: array of cell data

    CellData Format:
    float[]: psi(:, :, :, global cell index)
*/
void PsiData::writeToFile(const std::string &filename)
{
    MPI_File file;
    char outputName[32] = {
        'T', 'y', 'c', 'h', 'o', ' ', '2', ' ', 'P', 's', 'i', ' ',
        'O', 'u', 't', 'p', 'u', 't', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    uint64_t restOfHeader[4];


    // Fill in rest of header
    // version, number of cells, number of angles, and number of group
    restOfHeader[0] = 1;
    restOfHeader[1] = c_nc;
    Comm::gsum(restOfHeader[1]);
    restOfHeader[2] = c_na;
    restOfHeader[3] = c_ng;


    // Open file
    Comm::openFileForWrite(filename, file);


    // If rank 0, write header data
    if (Comm::rank() == 0) {
        double header[8];
        memcpy(header, outputName, 4 * sizeof(double));
        memcpy(&header[4], restOfHeader, 4 * sizeof(double));
        Comm::writeDoublesAt(file, 0, header, 8);
    }

   template <class T>
   std::vector <T>  dbldata(c_na*c_ng*c_nv);
    // Write data one cell at a time
    for (size_t cell = 0; cell < c_nc; cell++) {
        int dataSize = c_na * c_ng * c_nv;
        uint64_t globalCell = g_tychoMesh->getLGCell(cell);
        uint64_t offset = 8 + globalCell * dataSize;
        T *data = &c_data[index(0, 0, 0, cell)]; //not sure about this line
	         for(int i=0; i< dataSize; i++){
		           dbldata[i] = (double)data[i];
	}
        Comm::writeDoublesAt(file, offset, dbldata.data(), dataSize);
    }


    // Close file
    Comm::closeFile(file);
}
