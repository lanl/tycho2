
#include "SweeperPBJ.hh"
#include "Global.hh"
#include "TraverseGraph.hh"
#include "Priorities.hh"
#include "Transport.hh"
#include "PsiData.hh"
#include <omp.h>
#include <vector>
#include <math.h>


using namespace std;


/*
    populateLocalPsiBound
    
    Put data from neighboring cells into localPsiBound(fvrtx, face, group).
*/
static
void populateLocalPsiBound(const UINT angle, const UINT cell, 
                           const PsiData &psi, const PsiData &psiBound,
                           Mat3<double> &localPsiBound)
{
    // Default to 0.0
    localPsiBound = 0.0;
    
    // Populate if incoming flux
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_spTychoMesh->isIncoming(angle, cell, face)) {
            UINT neighborCell = g_spTychoMesh->getAdjCell(cell, face);
            
            // In local mesh
            if (neighborCell != TychoMesh::BOUNDARY_FACE) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT neighborVrtx = 
                        g_spTychoMesh->getNeighborVrtx(cell, face, fvrtx);
                    localPsiBound(fvrtx, face, group) = 
                        psi(neighborVrtx, angle, neighborCell, group);
                }
            }
            
            // Not in local mesh
            else if (g_spTychoMesh->getAdjRank(cell, face) != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT side = g_spTychoMesh->getSide(cell, face);
                    localPsiBound(fvrtx, face, group) = 
                        psiBound(side, angle, fvrtx, group);
                }
            }
        }
    }}
}



/*
    SweepData
    
    Holds psi and other data for the sweep.
*/
class SweepData : public TraverseData
{
public:
    
    /*
        Constructor
    */
    SweepData(PsiData &psi, const PsiData &source, PsiData &psiBound, 
              const double sigmaTotal)
    : c_psi(psi), c_psiBound(psiBound), 
      c_source(source), c_sigmaTotal(sigmaTotal)
    { }
    
    
    /*
        data
        
        Return face data for (cell, angle) pair.
    */
    virtual const char* getData(UINT cell, UINT face, UINT angle)
    {
        Mat2<double> localFaceData(g_nVrtxPerFace, g_nGroups);
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            UINT vrtx = g_spTychoMesh->getFaceToCellVrtx(cell, face, fvrtx);
            localFaceData(fvrtx, group) = c_psi(vrtx, angle, cell, group);
        }}
        
        return (char*) (&localFaceData[0]);
    }
    
    
    /*
        getDataSize
    */
    virtual size_t getDataSize()
    {
        return g_nGroups * g_nVrtxPerFace * sizeof(double);
    }
    
    
    /*
        sideData
        
        Set face data for (side, angle) pair.
    */
    virtual void setSideData(UINT side, UINT angle, const char *data)
    {
        Mat2<double> localFaceData(g_nVrtxPerFace, g_nGroups, (double*)data);
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            c_psiBound(side, angle, fvrtx, group) = localFaceData(fvrtx, group);
        }}
    }
    
    
    /*
        getPriority
        
        Return a priority for the cell/angle pair.
        Not needed for this class, so it is just set to a constant.
    */
    virtual UINT getPriority(UINT cell, UINT angle)
    {
        UNUSED_VARIABLE(cell);
        UNUSED_VARIABLE(angle);
        return 1;
    }
    
    
    /*
        update
        
        Updates psi for a given (cell, angle) pair.
    */
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell])
    {
        UNUSED_VARIABLE(adjCellsSides);
        UNUSED_VARIABLE(bdryType);
        
        Mat2<double> localSource(g_nVrtxPerCell, g_nGroups);
        Mat2<double> localPsi(g_nVrtxPerCell, g_nGroups);
        Mat3<double> localPsiBound(g_nVrtxPerFace, g_nFacePerCell, g_nGroups);
        
        
        // Populate localSource
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            localSource(vrtx, group) = c_source(vrtx, angle, cell, group);
        }}
        
        
        // Populate localPsiBound
        populateLocalPsiBound(angle, cell, c_psi, c_psiBound, localPsiBound);
        
        
        // Transport solve
        Transport::solve(cell, angle, c_sigmaTotal, 
                         localPsiBound, localSource, localPsi);
        
        
        // localPsi -> psi
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            c_psi(vrtx, angle, cell, group) = localPsi(vrtx, group);
        }}
    }
    
private:
    PsiData &c_psi;
    PsiData &c_psiBound;
    const PsiData &c_source;
    const double c_sigmaTotal;
};


/*
    MetaData struct
*/
struct MetaData
{
    UINT gSide;
    UINT angle;
    UINT cell;
    UINT face;
};


/*
    commSides
*/
void commSides(const vector<UINT> &adjRanks,
               const vector<vector<MetaData>> &sendMetaData,
               const vector<UINT> &numSendPackets,
               const vector<UINT> &numRecvPackets,
               SweepData &sweepData)
{
    int mpiError;
    UINT numToRecv;
    UINT numAdjRanks = adjRanks.size();
    vector<MPI_Request> mpiRecvRequests(numAdjRanks);
    vector<MPI_Request> mpiSendRequests(numAdjRanks);
    vector<vector<char>> dataToSend(numAdjRanks);
    vector<vector<char>> dataToRecv(numAdjRanks);
    
    
    // Data structures to send/recv packets
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        UINT packetSize = 2 * sizeof(UINT) + sweepData.getDataSize();
        dataToSend[rankIndex].resize(packetSize * numSendPackets[rankIndex]);
        dataToRecv[rankIndex].resize(packetSize * numRecvPackets[rankIndex]);
    }
    
    
    // Irecv data
    numToRecv = 0;
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        
        if (dataToRecv[rankIndex].size() > 0) {
            int tag = 0;
            int adjRank = adjRanks[rankIndex];
            MPI_Request request;
            mpiError = MPI_Irecv(dataToRecv[rankIndex].data(), 
                                 dataToRecv[rankIndex].size(), 
                                 MPI_BYTE, adjRank, tag, MPI_COMM_WORLD,
                                 &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiRecvRequests[rankIndex] = request;
            numToRecv++;
        }
        
        else {
            mpiRecvRequests[rankIndex] = MPI_REQUEST_NULL;
        }
    }
    
    
    // Update data to send and Isend it
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        
        if (dataToSend[rankIndex].size() > 0) {
            for (UINT metaDataIndex = 0; 
                 metaDataIndex < sendMetaData[rankIndex].size(); 
                 metaDataIndex++)
            {
                UINT gSide = sendMetaData[rankIndex][metaDataIndex].gSide;
                UINT angle = sendMetaData[rankIndex][metaDataIndex].angle;
                UINT cell  = sendMetaData[rankIndex][metaDataIndex].cell;
                UINT face  = sendMetaData[rankIndex][metaDataIndex].face;
                const char *data = sweepData.getData(cell, face, angle);
            
                char *ptr = &dataToSend[rankIndex][metaDataIndex];
                memcpy(ptr, &gSide, sizeof(UINT));
                ptr += sizeof(UINT);
                memcpy(ptr, &angle, sizeof(UINT));
                ptr += sizeof(UINT);
                memcpy(ptr, data, sweepData.getDataSize());
            }
            
            int tag = 0;
            int adjRank = adjRanks[rankIndex];
            MPI_Request request;
            mpiError = MPI_Isend(dataToSend[rankIndex].data(), 
                                 dataToSend[rankIndex].size(), 
                                 MPI_BYTE, adjRank, tag, MPI_COMM_WORLD, 
                                 &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiSendRequests[rankIndex] = request;
        }
        
        else {
            mpiSendRequests[rankIndex] = MPI_REQUEST_NULL;
        }
    }
    
    
    // Get data from Irecv
    for (UINT numWaits = 0; numWaits < numToRecv; numWaits++) {
        
        // Wait for a data packet to arrive
        int rankIndex;
        mpiError = MPI_Waitany(mpiRecvRequests.size(), mpiRecvRequests.data(), 
                               &rankIndex, MPI_STATUS_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
        
        
        // Process Data
        UINT packetSize = 2 * sizeof(UINT) + sweepData.getDataSize();
        UINT numPackets = dataToRecv[rankIndex].size() / packetSize;
        for (UINT packetIndex = 0; packetIndex < numPackets; packetIndex++) {
            char *ptr = &dataToRecv[rankIndex][packetIndex * packetSize];
            UINT gSide = 0;
            UINT angle = 0;
            memcpy(&gSide, ptr, sizeof(UINT));
            ptr += sizeof(UINT);
            memcpy(&angle, ptr, sizeof(UINT));
            ptr += sizeof(UINT);
            UINT side = g_spTychoMesh->getGLSide(gSide);
            sweepData.setSideData(side, angle, ptr);
        }
    }
    
    
    // Wait on send to complete
    if (mpiSendRequests.size() > 0) {
        mpiError = MPI_Waitall(mpiSendRequests.size(), mpiSendRequests.data(), 
                               MPI_STATUSES_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
    }
}


/*
    Constructor
*/
SweeperPBJ::SweeperPBJ(const double sigmaTotal)
{
    c_sigmaTotal = sigmaTotal;
}


/*
    sweep
*/
void SweeperPBJ::sweep(PsiData &psi, const PsiData &source)
{
    const bool doComm = false;
    const UINT maxComputePerStep = UINT_MAX;
    const UINT maxIter = 100;
    const double tolerance = 1e-5;
    
    
    // Set initial guess for psiBound
    PsiData psiBound(g_spTychoMesh->getNSides(), g_quadrature->getNumAngles(), 
                     g_nVrtxPerFace, g_nGroups);
    psiBound.setToValue(0.0);  // Change to something more reasonable.
    
    
    // Set psi0
    PsiData psi0(g_nVrtxPerCell, g_quadrature->getNumAngles(), 
                 g_spTychoMesh->getNCells(), g_nGroups);
    
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
        psi0(vertex, angle, cell, group) = psi(vertex, angle, cell, group);
    }}}}
    
    
    // Create SweepData for traversal
    SweepData sweepData(psi, source, psiBound, c_sigmaTotal);
    
    
    // Get adjacent ranks
    vector<UINT> adjRanks;
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        UINT adjRank = g_spTychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
        
        if (adjCell == TychoMesh::BOUNDARY_FACE && 
            adjRank != TychoMesh::BAD_RANK &&
            std::count(adjRanks.begin(), adjRanks.end(), adjRank) == 0)
        {
            adjRanks.push_back(adjRank);
        }
    }}
    
    
    // Populate sendMetaData, numSendPackets, and numRecvPackets
    vector<vector<MetaData>> sendMetaData(adjRanks.size());
    vector<UINT> numSendPackets(adjRanks.size());
    vector<UINT> numRecvPackets(adjRanks.size());
    
    for (UINT rankIndex = 0; rankIndex < adjRanks.size(); rankIndex++) {
        
        numSendPackets[rankIndex] = 0;
        numRecvPackets[rankIndex] = 0;
        
        for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
        for (UINT face = 0; face < g_nFacePerCell; face++) {
        
            UINT adjRank = g_spTychoMesh->getAdjRank(cell, face);        
            if (adjRank == adjRanks[rankIndex]) {
                for (UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
                    if (g_spTychoMesh->isOutgoing(angle, cell, face)) {
                        MetaData md;
                        UINT side = g_spTychoMesh->getSide(cell, face);
                        md.gSide = g_spTychoMesh->getLGSide(side);
                        md.angle = angle;
                        md.cell  = cell;
                        md.face  = face;
                        sendMetaData[rankIndex].push_back(md);
                        
                        numSendPackets[rankIndex]++;
                    }
                    else {
                        numRecvPackets[rankIndex]++;
                    }
                }
            }
        }}
    }
    
    
    // Sweep till converged
    UINT iter = 0;
    while (iter < maxIter) {
        
        // Sweep
        traverseGraph(maxComputePerStep, sweepData, doComm, MPI_COMM_WORLD, 
                      Direction_Forward);
        
        // Check tolerance and psi0 = psi1
        double errL1 = 0.0;
        double normL1 = 0.0;
        for (UINT group = 0; group < g_nGroups; ++group) {
        for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
        for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
            double p0 = psi0(vertex, angle, cell, group);
            double p1 = psi(vertex, angle, cell, group);
            
            errL1  += fabs(p0 - p1);
            normL1 += fabs(p1);
            
            psi0(vertex, angle, cell, group) = p1;
        }}}}
        
        if (errL1 / normL1 < tolerance)
            break;
        
        
        // Communicate
        commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, 
                  sweepData);
    }
}

