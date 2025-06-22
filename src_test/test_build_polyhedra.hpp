#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"
#include <string>
#include <cmath>
#include <list>
#include <map>
#include <vector>
#include <set>

#include <gtest/gtest.h>
#include "Import.hpp"
#include "Export.hpp"
#include "Triangulation.hpp"
#include "Dual.hpp"
#include "ShortPath.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

namespace PolyhedralLibrary{
	
	/*PolyhedralMesh OldMesh;

	OldMesh.NumCell0Ds = 6;
	OldMesh.Cell0DsId = {0,1,2,3,4,5};

	OldMesh.Cell0DsCoordinates = MatrixXd::Zero(3, OldMesh.NumCell0Ds);
	OldMesh.Cell0DsCoordinates.col(0) = Vector3d(1, 0, 0);
	OldMesh.Cell0DsCoordinates.col(1) = Vector3d(-1, 0, 0);
	OldMesh.Cell0DsCoordinates.col(2) = Vector3d(0, 1, 0);
	OldMesh.Cell0DsCoordinates.col(3) = Vector3d(0, -1, 0);
	OldMesh.Cell0DsCoordinates.col(4) = Vector3d(0, 0, 1);
	OldMesh.Cell0DsCoordinates.col(5) = Vector3d(0, 0, -1);

	OldMesh.ShortPathCell0Ds[0] = {0,1,2,3,4,5};

	OldMesh.NumCell1Ds = 12;
	OldMesh.Cell1DsId = {0,1,2,3,4,5,6,7,8,9,10,11};

	OldMesh.Cell1DsExtrema = MatrixXi::Zero(2, OldMesh.NumCell1Ds);
	OldMesh.Cell1DsExtrema.col(0) = Vector2i(0,2);
	OldMesh.Cell1DsExtrema.col(1) = Vector2i(0,3);
	OldMesh.Cell1DsExtrema.col(2) = Vector2i(0,4);
	OldMesh.Cell1DsExtrema.col(3) = Vector2i(0,5);
	OldMesh.Cell1DsExtrema.col(4) = Vector2i(1,2);
	OldMesh.Cell1DsExtrema.col(5) = Vector2i(1,3);
	OldMesh.Cell1DsExtrema.col(6) = Vector2i(1,4);
	OldMesh.Cell1DsExtrema.col(7) = Vector2i(1,5);
	OldMesh.Cell1DsExtrema.col(8) = Vector2i(2,4);
	OldMesh.Cell1DsExtrema.col(9) = Vector2i(2,5);
	OldMesh.Cell1DsExtrema.col(10) = Vector2i(3,4);
	OldMesh.Cell1DsExtrema.col(11) = Vector2i(3,5);

	OldMesh.ShortPathCell1Ds[0] = {0,1,2,3,4,5,6,7,8,9,10,11};

	OldMesh.NumCell2Ds = 8;
	OldMesh.Cell2DsId = {0,1,2,3,4,5,6,7};
	OldMesh.Cell2DsVertices = {{0,2,4},{2,1,4},{1,3,4},{3,0,4},{0,2,5},{2,1,5},{1,3,5},{3,0,5}};
	OldMesh.Cell2DsEdges = {{0,8,2},{4,6,8},{5,10,6},{1,2,10},{0,9,3},{4,7,9},{5,11,7},{1,3,11}};

	OldMesh.NumCell3Ds = 1;
	OldMesh.Cell3DsId = {0};
	OldMesh.Cell3DsVertices = {{0,1,2,3,4,5}};
	OldMesh.Cell3DsEdges = {{0,1,2,3,4,5,6,7,8,9,10,11}};
	OldMesh.Cell3DsFaces = {{0,1,2,3,4,5,6,7}};
	*/
	
	//con questa inizializzazione potrebbe funzionare, ma non la prende nei test...forse va inserita dentro ogni test? forse conviene usare import a questo punto... o comunque provare con figure 2D quando si può e non solidi veri e propri

	PolyhedralMesh NewMesh;
	
TEST(TestChecking, TestCheckAddEdges1)
{
	PolyhedralMesh mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	Vector2i edge = {0,5};
	int id_prova = CheckAddEdges(mesh,edge,mesh.NumCell1Ds);
	int id_corretto = 3; //3 è l'id corrispondente al lato con gli estremi indicati in id_prova
	EXPECT_EQ(id_prova,id_corretto);
}

/*TEST(TestChecking, TestCheckAddEdges1b)
{
    PolyhedralMesh mesh;
    string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
    string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
    string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
    string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
    ImportMesh(mesh, file0Ds, file1Ds, file2Ds, file3Ds);
    int originalNumEdges = mesh.NumCell1Ds;
    Vector2i edge = {4, 5}; //nuovo lato da aggiungere, non presente
    // Pre-allocazione se necessario (se CheckAddEdges assume spazio sufficiente)
    if (mesh.Cell1DsExtrema.cols() <= originalNumEdges) {
        mesh.Cell1DsExtrema.conservativeResize(2, originalNumEdges + 10);
    }
    int id_edge = originalNumEdges - 1; 
    int returnedId = CheckAddEdges(mesh, edge, id_edge);
    EXPECT_EQ(returnedId, originalNumEdges);
    EXPECT_EQ(mesh.Cell1DsId.back(), originalNumEdges);
    EXPECT_EQ(mesh.Cell1DsExtrema(0, originalNumEdges), 4);
    EXPECT_EQ(mesh.Cell1DsExtrema(1, originalNumEdges), 5);
}


TEST(TestChecking, TestCheckAddEdges2b)
{
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	Vector2i edge = {4,5}; //non è tra quelli già presenti nel file di input
	OldMesh.Cell1DsId.resize(1,OldMesh.NumCell0Ds + 1);
	OldMesh.Cell1DsExtrema.conservativeResize(2,OldMesh.NumCell0Ds + 1);
	int id_prova = CheckAddEdges(OldMesh, edge, OldMesh.NumCell1Ds);
	int id_corretto = OldMesh.NumCell1Ds + 1;
	EXPECT_EQ(id_prova,id_corretto);
}*/

//si potrebbero fare questi test sopra anche senza file di input, ma con una figura geometrica a caso, con vertici semplici

TEST(TestChecking, TestCheckAddEdges2)
{
	//aggiungo un nuovo lato che non esiste
    PolyhedralMesh mesh;
    mesh.Cell1DsExtrema = MatrixXi::Zero(2,10);
	mesh.Cell1DsExtrema.col(0) = Vector2i(0,1);
	mesh.Cell1DsExtrema.col(1) = Vector2i(1,2);
	mesh.Cell1DsExtrema.col(2) = Vector2i(2,0);
	mesh.Cell1DsId = {0,1,2};
	int id_iniziale = 2;
    Vector2i edge(2,3);
    int id_prova = CheckAddEdges(mesh, edge, id_iniziale);
    EXPECT_EQ(id_prova, 3);              
}


TEST(TestChecking, TestCheckAddVertices1)
{
	PolyhedralMesh mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	Vector3d vertex={0,0,1};
	int id_prova=CheckAddVertices(mesh,vertex,mesh.NumCell0Ds);
	int id_corretto=4; //4 è l'id corrispondente al vertice con gli estremi indicati in id_prova
	EXPECT_EQ(id_prova,id_corretto);
}

TEST(TestChecking, TestCheckAddVertices2)
{
	//aggiungo un nuovo vertice che non esiste
	PolyhedralMesh mesh;
    mesh.Cell0DsCoordinates = MatrixXd::Zero(3,10);
	mesh.Cell0DsCoordinates.col(0) = Vector3d(1, 0, 1);
	mesh.Cell0DsCoordinates.col(1) = Vector3d(-1, -1, 0);
	mesh.Cell0DsCoordinates.col(2) = Vector3d(0, 1, 1);
	mesh.Cell0DsId={0,1,2};
	int id_iniziale=2;
	Vector3d vertex{0,1,0};
	int id_prova=CheckAddVertices(mesh,vertex,id_iniziale);
	EXPECT_EQ(id_prova,3);
}

TEST(TestChecking, TestTriangulationTypeI)
{
	// confrontro tra numero di vertici, lati e facce del poliedro triangolato (ootenuto con la funzione) e i valori che si ottengono applicando le formule
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	/*OldMesh.Cell0DsCoordinates = MatrixXd::Zero(3,10);
	OldMesh.Cell0DsCoordinates.col(0) = Vector3d(0, 0, 0);
	OldMesh.Cell0DsCoordinates.col(1) = Vector3d(1, 0, 0);
	OldMesh.Cell0DsCoordinates.col(2) = Vector3d(0.5, 0.8660254038, 0);
	OldMesh.Cell2DsVertices = {{0,1,2}};*/
	PolyhedralMesh NewMesh;
	int p=3;
	int q=4;
	int n=2;
	TriangulationTypeI(OldMesh,NewMesh,p,q,n); 
	PolyhedralMesh mesh;
	mesh.Cell3DsId = {0};
	mesh.Cell3DsVertices = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17}};
	mesh.Cell3DsEdges = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47}};
	mesh.Cell3DsFaces = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}};
	//mesh.Cell3DsVertices = {{0,1,2,3,4,5,6,7,8,9}};
	//mesh.Cell3DsEdges = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17}};
	//mesh.Cell3DsFaces = {{0}};
	EXPECT_EQ(NewMesh.Cell3DsId,mesh.Cell3DsId);
	EXPECT_EQ(NewMesh.Cell3DsVertices,mesh.Cell3DsVertices);
	EXPECT_EQ(NewMesh.Cell3DsEdges,mesh.Cell3DsEdges);
	EXPECT_EQ(NewMesh.Cell3DsFaces,mesh.Cell3DsFaces);
}

TEST(TestChecking, TestTriangulationTypeII)
{
	// confrontro tra numero di vertici, lati e facce del poliedro triangolato (ootenuto con la funzione) e i valori che si ottengono applicando le formule
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	/*OldMesh.NumCell0Ds = 6;
	OldMesh.NumCell1Ds = 12;
	OldMesh.NumCell2Ds = 8;
	OldMesh.Cell0DsCoordinates = MatrixXd::Zero(3,10);
	OldMesh.Cell0DsCoordinates.col(0) = Vector3d(0, 0, 0);
	OldMesh.Cell0DsCoordinates.col(1) = Vector3d(1, 0, 0);
	OldMesh.Cell0DsCoordinates.col(2) = Vector3d(0.5, 0.8660254038, 0);
	OldMesh.Cell2DsVertices = {{0,1,2}};*/
	PolyhedralMesh NewMesh;
	//int p=3;
	//int q=4;
	int n=2;
	TriangulationTypeII(OldMesh,NewMesh,n);
	PolyhedralMesh mesh;
	mesh.Cell3DsId = {0};
	mesh.Cell3DsVertices = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73}};
	mesh.Cell3DsEdges = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215}};
	mesh.Cell3DsFaces = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143}};
	EXPECT_EQ(NewMesh.Cell3DsId,mesh.Cell3DsId);
	EXPECT_EQ(NewMesh.Cell3DsVertices,mesh.Cell3DsVertices);
	EXPECT_EQ(NewMesh.Cell3DsEdges,mesh.Cell3DsEdges);
	EXPECT_EQ(NewMesh.Cell3DsFaces,mesh.Cell3DsFaces);
}

TEST(TestChecking, TestDualConstructor)
{
	// confrontro tra il numero di vertici, lati e facce del poliedro duale che si ottine utilizzando la funzione e i valori che si ottengono applicando le formule
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh DualMesh;
	DualConstructor(OldMesh,DualMesh);
	PolyhedralMesh mesh;
	mesh.Cell3DsId = {0};
	mesh.Cell3DsVertices = {{0,1,2,3,4,5,6,7}};
	mesh.Cell3DsEdges = {{0,1,2,3,4,5,6,7,8,9,10,11}};
	mesh.Cell3DsFaces = {{0,1,2,3,4,5}};
	EXPECT_EQ(DualMesh.Cell3DsId,mesh.Cell3DsId);
	EXPECT_EQ(DualMesh.Cell3DsVertices,mesh.Cell3DsVertices);
	EXPECT_EQ(DualMesh.Cell3DsEdges,mesh.Cell3DsEdges);
	EXPECT_EQ(DualMesh.Cell3DsFaces,mesh.Cell3DsFaces);
}

TEST(TestChecking, TestFindEdge)
{
	PolyhedralMesh mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	int v0=2;
	int v1=4;
	int edgeId=FindEdge(mesh,v0,v1);
	EXPECT_EQ(edgeId,8);
}

TEST(TestChecking, TestCreateAdjacencyList)
{
	PolyhedralMesh mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	vector<list<int>> adjList(mesh.NumCell0Ds);
	CreateAdjacencyList(mesh,adjList);
	vector<list<int>> lista_corretta= {{2,3,4,5},{2,3,4,5},{0,1,4,5},{0,1,4,5},{0,1,2,3},{0,1,2,3}};
	EXPECT_EQ(adjList,lista_corretta);
}

TEST(TestChecking, TestCreateWeightsMatrix)
{
	PolyhedralMesh mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	MatrixXd weightsEdges = MatrixXd::Ones(mesh.NumCell0Ds,mesh.NumCell0Ds);
	CreateWeightsMatrix(mesh,weightsEdges);
	MatrixXd matrice_corretta = MatrixXd::Ones(mesh.NumCell0Ds,mesh.NumCell0Ds);
	matrice_corretta << 
		1.0,1.0,1.4142135623730951,1.4142135623730951,1.4142135623730951,1.4142135623730951,
		1.0,1.0,1.4142135623730951,1.4142135623730951,1.4142135623730951,1.4142135623730951,
		1.4142135623730951,1.4142135623730951,1.0,1.0,1.4142135623730951,1.4142135623730951,
		1.4142135623730951,1.4142135623730951,1.0,1.0,1.4142135623730951,1.4142135623730951,
		1.4142135623730951,1.4142135623730951,1.4142135623730951,1.4142135623730951,1.0,1.0,
		1.4142135623730951,1.4142135623730951,1.4142135623730951,1.4142135623730951,1.0,1.0;
	EXPECT_TRUE(weightsEdges.isApprox(matrice_corretta, 1e-15));
}

TEST(TestChecking, TestComputeDistances)
{
	PolyhedralMesh mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	vector<list<int>> adjacencyList(mesh.NumCell0Ds);
	CreateAdjacencyList(mesh, adjacencyList);
	MatrixXd weightsEdges = MatrixXd::Ones(mesh.NumCell0Ds,mesh.NumCell0Ds)*INFINITY;
	CreateWeightsMatrix(mesh, weightsEdges);
	vector<int> predecessors(mesh.NumCell0Ds);
	vector<double> distances(mesh.NumCell0Ds);
	int sourceNode=3;
	int destinationNode=15;
	ComputeDistances(adjacencyList, sourceNode, destinationNode, weightsEdges, mesh.NumCell0Ds, predecessors, distances);
	double distanza_corretta=3.21143e-322;
	EXPECT_EQ(distances[destinationNode],distanza_corretta);
}

/*TEST(TestChecking, TestFindShortestPath)
{
	// controllo del cammino minimo
	PolyhedralMesh Mesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(Mesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh NewMesh;
	int p=3;
	int q=4;
	int n=2;
	TriangulationTypeI(Mesh,NewMesh,p,q,n);
	int nodo_partenza=3;
	int nodo_arrivo=15;
	unsigned int numero_lati=0;
	double lunghezza_cammino=0.0;
	FindShortestPath(NewMesh,nodo_partenza,nodo_arrivo,numero_lati,lunghezza_cammino);
	PolyhedralMesh mesh;
	mesh.ShortPathCell0Ds[1] = {3,1,2,5,15};
	mesh.ShortPathCell1Ds[1] = {3,1,8,34};
	unsigned int num_lati_corretto=5;
	double lung_cammino_corretta=26.23;
	EXPECT_EQ(NewMesh.ShortPathCell0Ds[1],mesh.ShortPathCell0Ds[1]);
	EXPECT_EQ(NewMesh.ShortPathCell1Ds[1],mesh.ShortPathCell1Ds[1]);
	EXPECT_EQ(numero_lati,num_lati_corretto);
	EXPECT_EQ(lunghezza_cammino,lung_cammino_corretta);
}*/
}