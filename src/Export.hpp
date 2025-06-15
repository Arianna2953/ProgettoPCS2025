#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
// Export the Cell0D properties from PolyhedralMesh
// mesh: a PolyhedralMesh struct
// return the result of the creation of the new polyhedron, true if is success, false otherwise
bool Exportfile0Ds(const PolyhedralMesh & polyNew);

// Export the Cell1D properties from PolyhedralMesh
// mesh: a PolyhedralMesh struct
// return the result of the creation of the new polyhedron, true if is success, false otherwise
bool Exportfile1Ds(const PolyhedralMesh & polyNew);

// Export the Cell2D properties from PolyhedralMesh
// mesh: a PolyhedralMesh struct
// return the result of the creation of the new polyhedron, true if is success, false otherwise
bool Exportfile2Ds(const PolyhedralMesh & polyNew);

// Export the Cell3D properties from PolyhedralMesh
// mesh: a PolyhedralMesh struct
// return the result of the creation of the new polyhedron, true if is success, false otherwise
bool Exportfile3Ds(const PolyhedralMesh & polyNew);

// General function for exporting the polyhedron from PolyhedralMesh
bool ExportPolyhedron(const PolyhedralMesh & polyNew);

}
