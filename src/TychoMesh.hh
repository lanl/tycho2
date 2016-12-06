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

/*   
    Notation convention for the mesh   

    Cell: global index for a Tetrahedron
    Node: global index of a Vertex
    Side: global index of a face side
          each face possibly has two sides
          if a face is on a boundary, then it contains a side
    Vrtx: locally defined with respect to either a cell or face
          if cell: index 0,1,2,3   if face: index 0,1,2
    Face: local indexing of a face from a cell: 0,1,2,3
*/

#ifndef __TYCHO_MESH_HH__
#define __TYCHO_MESH_HH__

#include "Mat.hh"
#include "Global.hh"
#include "Assert.hh"
#include <map>


class TychoMesh 
{
public:
    // Constructor
    TychoMesh(const std::string &filename);
    
    
    // Get data
    UINT getNSides() const
        { return c_nSides; }
    UINT getNNodes() const
        { return c_nNodes; }
    double getNodeCoord(UINT node, UINT dim) const
        { return c_nodeCoords(node, dim); }
    UINT getCellNode(UINT cell, UINT vrtx) const
        { return c_cellNodes(cell, vrtx); }
    UINT getAdjCell(UINT cell, UINT face) const
        { return c_adjCell(cell, face); }
    UINT getAdjFace(UINT cell, UINT face) const
        { return c_adjCell(cell, face); }
    UINT getSideCell(UINT side) const
        { return c_sideCell(side); }
    UINT getSide(UINT cell, UINT face) const
        { return c_side(cell, face); }
    UINT getLGSide(const UINT side) const
        { return c_lGSides(side); }
    UINT getGLSide(const UINT side) const
        { return c_gLSides.find(side)->second; }
    UINT getLGCell(const UINT cell) const
        { return c_lGCells(cell); }
    UINT getAdjRank(const UINT cell, const UINT face) const
        { return c_adjProc(cell, face); }
    UINT getFaceToCellVrtx(const UINT cell, const UINT face, const UINT fvrtx) const
        { return c_faceToCellVrtx(cell, face, fvrtx); }
    UINT getCellToFaceVrtx(const UINT cell, const UINT face, const UINT cvrtx) const
        { Assert(cvrtx != face);
          return c_cellToFaceVrtx(cell, face, cvrtx); }
    double getOmegaDotN(UINT angle, UINT cell, UINT face) const
        { return c_omegaDotN(angle, cell, face); }
    double getCellVolume(const UINT cell) const
        { return c_cellVolume(cell); }
    double getFaceArea(const UINT cell, const UINT face) const
        { return c_faceArea(cell, face); }
    UINT getNeighborVrtx(const UINT cell, const UINT face, const UINT fvrtx) const
        { return c_neighborVrtx(cell, face, fvrtx); }
    bool isOutgoing(const UINT angle, const UINT cell, const UINT face) const
        { return getOmegaDotN(angle, cell, face) > 0; }
    bool isIncoming(const UINT angle, const UINT cell, const UINT face) const
        { return !isOutgoing(angle, cell, face); }
    UINT getAdjCellFromSide(const UINT side) const
        { return c_adjCellFromSide(side); }
    UINT getAdjFaceFromSide(const UINT side) const
        { return c_adjFaceFromSide(side); }
    UINT getCellMaterial(const UINT cell) const
        { return c_cellMaterial(cell); }
    
    
    // Arbitrary value to mark any face that lies on a boundary.
    static const UINT BOUNDARY_FACE = UINT64_MAX;
    static const UINT NOT_BOUNDARY_FACE = UINT64_MAX;
    static const UINT BAD_RANK = UINT64_MAX;
    

    // Structures
    struct FaceCoords
    {
        double c[3][3];
    };

    struct CellCoords
    {
        double c[4][3];
    };
    
private:
    void readTychoMesh(const std::string &filename);
    CellCoords getCellVrtxCoords(UINT cell) const;
    FaceCoords getFaceVrtxCoords(UINT cell, UINT face) const;
    UINT getCellVrtx(const UINT cell, const UINT node) const;
    
    UINT c_nSides;
    UINT c_nNodes;
    Mat2<double> c_nodeCoords;      // (node, dim) -> coord
    Mat2<UINT> c_cellNodes;         // (cell, vrtx) -> node
    Mat2<UINT> c_adjCell;           // (cell, face) -> cell
    Mat2<UINT> c_adjFace;           // (cell, face) -> face
    Mat1<UINT> c_sideCell;          // side -> cell
    Mat2<UINT> c_side;              // (cell, face) -> side
    Mat1<UINT> c_lGSides;           // local to global side numbering.
    Mat1<UINT> c_lGCells;           // local to global side numbering.
    std::map<UINT, UINT> c_gLSides; // global to local side numbering.
    Mat2<UINT> c_adjProc;           // (cell, face) -> adjacent proc
    Mat3<double> c_omegaDotN;       // (angle, cell, face) -> omega dot n
    Mat1<double> c_cellVolume;      // cell -> volume
    Mat2<double> c_faceArea;        // (cell, face) -> area
    Mat3<UINT> c_faceToCellVrtx;    // (cell, face, fvrtx) -> cvrtx
    Mat3<UINT> c_cellToFaceVrtx;    // (cell, face, cvrtx) -> fvrtx
    Mat3<UINT> c_neighborVrtx;      // (cell, face, fvrtx) -> vrtx
    Mat1<UINT> c_adjCellFromSide;   // side -> adj cell
    Mat1<UINT> c_adjFaceFromSide;   // side -> adj face
    Mat1<UINT> c_cellMaterial;      // cell -> material index
};


#endif
