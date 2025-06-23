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
    int id_prova = CheckAddEdges(mesh,edge,id_iniziale);
    EXPECT_EQ(id_prova,3);              
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
	//confrontro tra numero di vertici, lati e facce del poliedro triangolato (ottenuto con la funzione) e i valori che si ottengono applicando le formule
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
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
	EXPECT_EQ(NewMesh.Cell3DsId,mesh.Cell3DsId);
	EXPECT_EQ(NewMesh.Cell3DsVertices,mesh.Cell3DsVertices);
	EXPECT_EQ(NewMesh.Cell3DsEdges,mesh.Cell3DsEdges);
	EXPECT_EQ(NewMesh.Cell3DsFaces,mesh.Cell3DsFaces);
}

TEST(TestChecking, TestTriangulationTypeII)
{
	//confrontro tra numero di vertici, lati e facce del poliedro triangolato (ottenuto con la funzione) e i valori che si ottengono applicando le formule
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh NewMesh;
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
	//confrontro tra il numero di vertici, lati e facce del poliedro duale che si ottine utilizzando la funzione e i valori che si ottengono applicando le formule
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
	//verifico che il lato che unisce due vertici, trovato dalla funzione, corrisponda con il lato che li unisce realmente
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
	//verifico che la lista di adiacenza, trovata con la funzione, sia la medesima di quella trovata manualmente
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
	//verifico che la matrice dei pesi, che si crea tramite alla funzione, si differenzia dalla matrice calcolata manualmente di una tolleranza minima
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
	//verifico che la distanza tra i nodi, trovata con la funzione, sia la medesima di quella trovata manualmente
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
	ComputeDistances(adjacencyList, sourceNode, destinationNode, weightsEdges, predecessors, distances);
	double distanza_corretta=3.21143e-322;
	EXPECT_EQ(distances[destinationNode],distanza_corretta);
}

bool isApprox(double a, double b, double tolerance = 1e-15) {
  return abs(a - b) <= tolerance;
}

TEST(TestChecking, TestFindShortestPath)
{
	//verifico che vertici, lati e il numero di questi ultimo, trovati dalla funzione, corrispondono con quelli trovati manualmente; inoltre verifico anche che la lunghezza del cammino minimo, trovata dalla funzione, si differenzia dalla lunghezza calcolata manualmente di una tolleranza minima
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
	mesh.ShortPathCell0Ds[1] = {3,4,15};
	mesh.ShortPathCell1Ds[1] = {4,33};
	unsigned int num_lati_corretto=2;
	double lung_cammino_corretta= 1.7653668647301792e+00;
	EXPECT_EQ(NewMesh.ShortPathCell0Ds[1],mesh.ShortPathCell0Ds[1]);
	EXPECT_EQ(NewMesh.ShortPathCell1Ds[1],mesh.ShortPathCell1Ds[1]);
	EXPECT_EQ(numero_lati,num_lati_corretto);
	EXPECT_TRUE(isApprox(lunghezza_cammino,lung_cammino_corretta, 1e-15));
}
}