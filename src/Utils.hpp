#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
//Triangulate the faces of the polyhedron
//mesh: a PolyhedralMesh struct
//b & c: input numbers 
void TriangulateFaces(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& p,const int& q, const int& b,const int& c);

// Projects the vertices of the triangulation onto a unit sphere centered at the origin
// vertices: a vector containing 3 coordinates of the vertex
// return the normalized vector
bool projectOntoUnitSphere(Vector3d& v);

void TriangulationTypeI(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& p, const int& q, const int& n); 

int CheckAddEdges(PolyhedralMesh& poly, const Vector2i edge, int& id_edge);
// finire e correggere i commenti!!!!

}