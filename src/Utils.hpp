#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
// Import the polyhedral mesh and test if the mesh is correct
// mesh: a PolyhedralMesh struct
// return the result of the reading, true if is success, false otherwise
bool ImportMesh(PolyhedralMesh& mesh,const string& file0Ds,const string& file1Ds,const string& file2Ds,const string& file3Ds);

// Import the Cell0D properties from Cell0Ds_type.csv file
// mesh: a PolyhedralMesh struct
// return the result of the reading, true if is success, false otherwise
bool ImportCell0Ds(const string& file0Ds, PolyhedralMesh& mesh);

// Import the Cell1D properties from Cell1Ds_type.csv file
// mesh: a PolyhedralMesh struct
// return the result of the reading, true if is success, false otherwise
bool ImportCell1Ds(const string& file1Ds, PolyhedralMesh& mesh);

// Import the Cell2D properties from Cell2Ds_type.csv file
// mesh: a PolyhdralMesh struct
// return the result of the reading, true if is success, false otherwise
bool ImportCell2Ds(const string& file2Ds, PolyhedralMesh& mesh);

// Import the Cell3D properties from Cell3Ds_type.csv file
// mesh: a PolyhdralMesh struct
// return the result of the reading, true if is success, false otherwise
bool ImportCell3Ds(const string& file3Ds, PolyhedralMesh& mesh);
}
