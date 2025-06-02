#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {

struct PolyhedralMesh
{
    int NumCell0Ds = 0; // numero di Cell0D
    vector<int> Cell0DsId = {}; // id vertici, (1 x NumberCell0D)
    MatrixXd Cell0DsCoordinates = {}; // coordinate dei punti, (3 x NumberCell0D (x,y))
    map<int, list<int>> ShortPathCell0Ds = {}; // Cell0D markers

    int NumCell1Ds = 0; // numero di Cell1D
    vector<int> Cell1DsId = {}; // id lati, (1 x NumberCell1D)
    MatrixXi Cell1DsExtrema = {}; // estremi dei lati, (2 x NumberCell1D) (fromId,toId)
    map<int, list<int>> ShortPathCell1Ds = {}; // Cell1D markers

    int NumCell2Ds = 0; // numero di Cell2D
    vector<int> Cell2DsId = {}; // id poligoni, (1 x NumberCell2D)
    vector<vector<int>> Cell2DsVertices = {}; // id vertici del poligono, (1 x NumberCell2DVertices[NumberCell2D])
    vector<vector<int>> Cell2DsEdges = {}; // id lati del poligono, (1 x NumberCell2DEdges[NumberCell2D])
	//map<int, list<int>> MarkerCell2Ds = {}; // Cell2D markers
	
	int NumCell3Ds = 0; // numero di Cell3D
	vector<int> Cell3DsId = {}; // id poligoni, (1 x NumberCell3D)
	vector<vector<int>> Cell3DsVertices = {}; // id vertici del poligono, (1 x NumberCell3DVertices[NumberCell3D])
    vector<vector<int>> Cell3DsEdges = {}; // id lati del poligono, (1 x NumberCell3DEdges[NumberCell3D])
	vector<vector<int>> Cell3DsFaces = {}; // id facce del poligono, (1 x NumberCell3Faces[NumberCell3D])
	//map<int, list<int>> MarkerCell3Ds = {}; // Cell3D markers
};
}