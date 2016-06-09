/*!
    \file   Comm.cc
    \author Kris Garrett
    \date   2016
    \brief  Implements several MPI requests.
 */

#include "Comm.hh"
#include "Assert.hh"
#include <mpi.h>
#include <omp.h>


static const int INT_TAG = 1;
static const int DOUBLE_TAG = 2;


namespace Comm
{

/*! 
    \brief Returns MPI rank.
*/
int rank()
{
    int rank = 0;
    int result = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Insist(result == MPI_SUCCESS, "Comm::rank MPI error.\n");
    return rank;
}


/*! 
    \brief Returns MPI rank.
*/
int rank(MPI_Comm comm)
{
    int rank = 0;
    int result = MPI_Comm_rank(comm, &rank);
    Insist(result == MPI_SUCCESS, "Comm::rank MPI error.\n");
    return rank;
}


/*!
    \brief Returns number of MPI ranks.
*/
int numRanks()
{
    int numRanks = 0;
    int result = MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    Insist(result == MPI_SUCCESS, "Comm::numRanks MPI error.\n");
    return numRanks;
}


/*!
    \brief Sums x from all ranks.
*/
void gsum(double &x)
{
    double send = x;
    double recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::gsum(double) MPI error.\n");
    x = recv;
}


/*!
    \brief Sums x from all ranks.
*/
void gsum(int &x)
{
    int send = x;
    int recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::gsum(int) MPI error.\n");
    x = recv;
}


/*!
    \brief Sums x from all ranks.
*/
void gsum(UINT &x)
{
    UINT send = x;
    UINT recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::gsum(UINT) MPI error.\n");
    x = recv;
}


/*!
    \brief Calculates max of x from all ranks.
*/
void gmax(double &x)
{
    double send = x;
    double recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::gmax(double) MPI error.\n");
    x = recv;
}


/*!
    \brief Calculates max of x from all ranks.
*/
void gmax(double &x, MPI_Comm comm)
{
    double send = x;
    double recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, MPI_MAX, comm);
    Insist(result == MPI_SUCCESS, "Comm::gmax(double) MPI error.\n");
    x = recv;
}


/*!
    \brief Calculates max of x from all ranks.
*/
void gmax(int &x)
{
    int send = x;
    int recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::gmax(int) MPI error.\n");
    x = recv;
}


/*!
    \brief Calculates max of x from all ranks.
*/
void gmax(UINT &x)
{
    UINT send = x;
    UINT recv = 0;
    int result = MPI_Allreduce(&send, &recv, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::gmax(int) MPI error.\n");
    x = recv;
}


/*!
    \brief Blocking send of a single integer value.
*/
void sendInt(int i, int destination)
{
    int result = MPI_Send(&i, 1, MPI_INT, destination, INT_TAG, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::sendInt MPI error.\n");
}


/*!
    \brief Blocking send of a single integer value.
*/
void sendUInt(UINT i, int destination)
{
    int result = MPI_Send(&i, 1, MPI_UINT64_T, destination, INT_TAG, MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::sendUInt MPI error.\n");
}


/*!
    \brief Blocking send of a vector of integer values.
*/
void sendIntVector(const std::vector<int> &buffer, int destination)
{
    int *buffer1 = const_cast<int*>(buffer.data());
    int result = MPI_Send(buffer1, buffer.size(), MPI_INT, destination, INT_TAG, 
                          MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::sendIntVector MPI error.\n");
}


/*!
    \brief Blocking send of a vector of integer values.
*/
void sendUIntVector(const std::vector<UINT> &buffer, int destination)
{
    UINT *buffer1 = const_cast<UINT*>(buffer.data());
    int result = MPI_Send(buffer1, buffer.size(), MPI_UINT64_T, destination, INT_TAG, 
                          MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::sendUIntVector MPI error.\n");
}


/*!
    \brief Asynchronous send of int vector for a given tag.
*/
void iSendIntVector(const std::vector<int> &buffer, int destination, int tag, 
                    MPI_Request &request)
{
    int *buffer1 = const_cast<int*>(buffer.data());
    int result = MPI_Isend(buffer1, buffer.size(), MPI_INT, destination, tag, 
                           MPI_COMM_WORLD, &request);
    Insist(result == MPI_SUCCESS, "Comm::iSendIntVector MPI error.\n");
}


/*!
    \brief Asynchronous send of int vector for a given tag.
*/
void iSendUIntVector(const std::vector<UINT> &buffer, int destination, int tag, 
                    MPI_Request &request)
{
    UINT *buffer1 = const_cast<UINT*>(buffer.data());
    int result = MPI_Isend(buffer1, buffer.size(), MPI_UINT64_T, destination, tag, 
                           MPI_COMM_WORLD, &request);
    Insist(result == MPI_SUCCESS, "Comm::iSendUIntVector MPI error.\n");
}


/*!
    \brief Asynchronous send of double vector for a given tag.
*/
void iSendDoubleVector(const std::vector<double> &buffer, int destination, int tag, 
                       MPI_Request &request)
{
    double *buffer1 = const_cast<double*>(buffer.data());
    int result = MPI_Isend(buffer1, buffer.size(), MPI_DOUBLE, destination, tag, 
                           MPI_COMM_WORLD, &request);
    Insist(result == MPI_SUCCESS, "Comm::iSendDoubleVector MPI error.\n");
}


/*!
    \brief Blocking receive of a single integer value or a vector of integer values.
*/
void recvInt(int &i, int destination)
{
    int result = MPI_Recv(&i, 1, MPI_INT, destination, INT_TAG, MPI_COMM_WORLD, 
                          MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvInt MPI error.\n");
}


/*!
    \brief Blocking receive of a single integer value or a vector of integer values.
*/
void recvUInt(UINT &i, int destination)
{
    int result = MPI_Recv(&i, 1, MPI_UINT64_T, destination, INT_TAG, MPI_COMM_WORLD, 
                          MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvUInt MPI error.\n");
}


/*!
    \brief Blocking receive of a single integer value or a vector of integer values.
*/
void recvIntVector(std::vector<int> &buffer, int destination)
{
    int result = MPI_Recv(&buffer[0], buffer.size(), MPI_INT, destination, INT_TAG, 
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvIntVector MPI error.\n");
}


/*!
    \brief Blocking receive of a single integer value or a vector of integer values.
*/
void recvUIntVector(std::vector<UINT> &buffer, int destination)
{
    int result = MPI_Recv(&buffer[0], buffer.size(), MPI_UINT64_T, destination, INT_TAG, 
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvUIntVector MPI error.\n");
}


/*!
    \brief Blocking receive of int/double vector for a given tag.
*/
void recvIntVector(std::vector<int> &buffer, int destination, int tag)
{
    int result = MPI_Recv(&buffer[0], buffer.size(), MPI_INT, destination, tag, 
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvIntVector with tag MPI error.\n");
}


/*!
    \brief Blocking receive of int/double vector for a given tag.
*/
void recvUIntVector(std::vector<UINT> &buffer, int destination, int tag)
{
    int result = MPI_Recv(&buffer[0], buffer.size(), MPI_UINT64_T, destination, tag, 
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvUIntVector with tag MPI error.\n");
}


/*!
    \brief Blocking receive of int/double vector for a given tag.
*/
void recvDoubleVector(std::vector<double> &buffer, int destination, int tag)
{
    int result = MPI_Recv(&buffer[0], buffer.size(), MPI_DOUBLE, destination, tag, 
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::recvDoubleVector MPI error.\n");
}


/*!
    \brief Global MPI barrier.
*/
void barrier()
{
    int result = MPI_Barrier(MPI_COMM_WORLD);
    Insist(result == MPI_SUCCESS, "Comm::barrier MPI error.\n");
}


/*!
    \brief Open file for MPI-IO.
*/
void openFileForRead(const std::string &filename, MPI_File &file)
{
    char *filename1 = const_cast<char*>(filename.c_str());
    int result = MPI_File_open(MPI_COMM_WORLD, filename1, MPI_MODE_RDONLY, 
                               MPI_INFO_NULL, &file);
    Insist(result == MPI_SUCCESS, "Comm::openFileForRead MPI error.\n");
}


/*!
    \brief Close file for MPI-IO.
*/
void closeFile(MPI_File &file)
{
    int result = MPI_File_close(&file);
    Insist(result == MPI_SUCCESS, "Comm::closeFile MPI error.\n");
}


/*!
    \brief Seek to file position.
*/
void seek(const MPI_File &file, uint64_t position)
{
    MPI_Offset offset = position;
    int result = MPI_File_seek(file, offset, MPI_SEEK_SET);
    Insist(result == MPI_SUCCESS, "Comm::seek MPI error.\n");
}


/*!
    \brief Read a uint64_t from file.
*/
void readUint64(const MPI_File &file, uint64_t &data)
{
    int result = MPI_File_read(file, &data, sizeof(uint64_t), MPI_BYTE, 
                               MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::readSizeT MPI error.\n");
}


/*!
    \brief Read an array of uint64_t elements.
*/
void readUint64(const MPI_File &file, uint64_t *data, int numData)
{
    int result = MPI_File_read(file, data, sizeof(uint64_t) * numData, MPI_BYTE, 
                               MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::readSizeT[] MPI error.\n");
}


/*!
    \brief Read an array of double elements.
*/
void readDouble(const MPI_File &file, double *data, int numData)
{
    int result = MPI_File_read(file, data, sizeof(double) * numData, MPI_BYTE, 
                               MPI_STATUS_IGNORE);
    Insist(result == MPI_SUCCESS, "Comm::readDouble[] MPI error.\n");
}


} // End namespace Comm
