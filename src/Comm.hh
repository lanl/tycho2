/*
    Comm.hh
*/

#ifndef __COMM_HH__
#define __COMM_HH__

#include "Typedef.hh"
#include <mpi.h>
#include <vector>

namespace Comm
{

int rank();
int rank(MPI_Comm comm);
int numRanks();

void gsum(double &x);
void gsum(int &x);
void gsum(UINT &x);
void gmax(double &x);
void gmax(double &x, MPI_Comm comm);
void gmax(int &x);
void gmax(UINT &x);

void sendInt(int i, int destination);
void sendUInt(UINT i, int destination);
void sendIntVector(const std::vector<int> &buffer, int destination);
void sendUIntVector(const std::vector<UINT> &buffer, int destination);
void iSendIntVector(const std::vector<int> &buffer, int destination, int tag, 
                    MPI_Request &request);
void iSendUIntVector(const std::vector<UINT> &buffer, int destination, int tag, 
                     MPI_Request &request);
void iSendDoubleVector(const std::vector<double> &buffer, int destination, int tag, 
                       MPI_Request &request);

void recvInt(int &i, int destination);
void recvUInt(UINT &i, int destination);
void recvIntVector(std::vector<int> &buffer, int destination);
void recvUIntVector(std::vector<UINT> &buffer, int destination);
void recvIntVector(std::vector<int> &buffer, int destination, int tag);
void recvUIntVector(std::vector<UINT> &buffer, int destination, int tag);
void recvDoubleVector(std::vector<double> &buffer, int destination, int tag);

void barrier();

void openFileForRead(const std::string &filename, MPI_File &file);
void closeFile(MPI_File &file);
void seek(const MPI_File &file, uint64_t position);
void readUint64(const MPI_File &file, uint64_t &data);
void readUint64(const MPI_File &file, uint64_t *data, int numData);
void readDouble(const MPI_File &file, double *data, int numData);

} // End namespace

#endif
