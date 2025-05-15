#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
//Triangulate the faces of the polyhedron
//mesh: a PolyhedralMesh struct
//b & c: input numbers 
bool TriangulateFaces(PolyhedralMesh& mesh,const int& b,const int& c);

// Projects the vertices of the triangulation onto a unit sphere centered at the origin
// vertices: a vector containing 3 coordinates of the vertex
// return the normalized vector
void projectOntoUnitSphere(Vector3d& v);
}