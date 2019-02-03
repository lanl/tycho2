struct Tetrahedral {};

template<typename Args...>
class Type<Tetrahedral>{};


/*
constexpr static auto area(Face<Tetrahedral> face)
{
    double a[3], b[3];
    double crossProduct[3];
    
    a[0] = points.c[1][0] - points.c[0][0];
    a[1] = points.c[1][1] - points.c[0][1];
    a[2] = points.c[1][2] - points.c[0][2];
    
    b[0] = points.c[2][0] - points.c[0][0];
    b[1] = points.c[2][1] - points.c[0][1];
    b[2] = points.c[2][2] - points.c[0][2];
    
    auto c = cross_product(a, b);
    
    return 0.5 * sqrt(crossProduct[0] * crossProduct[0] + 
                      crossProduct[1] * crossProduct[1] + 
                      crossProduct[2] * crossProduct[2]);
}


static auto volume(Type<Tetrahedral> tet)
{
    double a[3], b[3], c[3];
    double crossProduct[3];
    
    a[0] = points.c[1][0] - points.c[0][0];
    a[1] = points.c[1][1] - points.c[0][1];
    a[2] = points.c[1][2] - points.c[0][2];
    
    b[0] = points.c[2][0] - points.c[0][0];
    b[1] = points.c[2][1] - points.c[0][1];
    b[2] = points.c[2][2] - points.c[0][2];
    
    c[0] = points.c[3][0] - points.c[0][0];
    c[1] = points.c[3][1] - points.c[0][1];
    c[2] = points.c[3][2] - points.c[0][2];
    
    getCrossProduct(a, b, crossProduct);
    
    return 1.0 / 6.0 * fabs(crossProduct[0] * c[0] + crossProduct[1] * c[1] + 
                            crossProduct[2] * c[2]);
}


static
vector<double> getNormal(TychoMesh::FaceCoords facePoints, 
                         TychoMesh::CellCoords cellPoints)
{
    double a[3], b[3];
    double crossProduct[3];
    double norm;
    UINT vertIndex;
    vector<double> normal(3);
    double dotProduct;
    
    
    // Get two vectors to cross
    a[0] = facePoints.c[1][0] - facePoints.c[0][0];
    a[1] = facePoints.c[1][1] - facePoints.c[0][1];
    a[2] = facePoints.c[1][2] - facePoints.c[0][2];
    
    b[0] = facePoints.c[2][0] - facePoints.c[0][0];
    b[1] = facePoints.c[2][1] - facePoints.c[0][1];
    b[2] = facePoints.c[2][2] - facePoints.c[0][2];
    
    
    // Get a normal vector
    getCrossProduct(a, b, crossProduct);
    
    norm = sqrt(crossProduct[0] * crossProduct[0] + 
                crossProduct[1] * crossProduct[1] + 
                crossProduct[2] * crossProduct[2]);
    
    normal[0] = crossProduct[0] / norm;
    normal[1] = crossProduct[1] / norm;
    normal[2] = crossProduct[2] / norm;
    
    
    // Find index of point in cellPoints not in facePoints
    for(vertIndex = 0; vertIndex < 4; vertIndex++) {
        bool pointEqual = false;
        for(UINT vert = 0; vert < 3; vert++) {
            if(cellPoints.c[vertIndex][0] == facePoints.c[vert][0] && 
               cellPoints.c[vertIndex][1] == facePoints.c[vert][1] && 
               cellPoints.c[vertIndex][2] == facePoints.c[vert][2])
            {
                pointEqual = true;
            }
        }
        
        if(pointEqual == false)
            break;
    }
    
    
    // Make sure normal vector is outward normal
    dotProduct = normal[0] * (cellPoints.c[vertIndex][0] - facePoints.c[0][0]) + 
                 normal[1] * (cellPoints.c[vertIndex][1] - facePoints.c[0][1]) + 
                 normal[2] * (cellPoints.c[vertIndex][2] - facePoints.c[0][2]);
    if(dotProduct > 0) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
    }
    
    
    return normal;
}
*/
