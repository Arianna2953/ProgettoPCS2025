#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "ImportExport.hpp"
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
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	Vector2i edge = {0,5};
	int id_prova = CheckAddEdges(OldMesh,edge,OldMesh.NumCell1Ds);
	int id_corretto = 3; //3 è l'id corrispondente al lato con gli estremi indicati in id_prova
	EXPECT_EQ(id_prova,id_corretto);
}

TEST(TestChecking, TestCheckAddEdges2)
{
    PolyhedralMesh OldMesh;
    string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
    string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
    string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
    string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";

    ImportMesh(OldMesh, file0Ds, file1Ds, file2Ds, file3Ds);

    int originalNumEdges = OldMesh.NumCell1Ds;
    Vector2i edge = {4, 5}; //nuovo lato da aggiungere, non presente

    // Pre-allocazione se necessario (se CheckAddEdges assume spazio sufficiente)
    if (OldMesh.Cell1DsExtrema.cols() <= originalNumEdges) {
        OldMesh.Cell1DsExtrema.conservativeResize(2, originalNumEdges + 10);
    }

    int id_edge = originalNumEdges - 1; 

    int returnedId = CheckAddEdges(OldMesh, edge, id_edge);

    EXPECT_EQ(returnedId, originalNumEdges);
    EXPECT_EQ(OldMesh.Cell1DsId.back(), originalNumEdges);
    EXPECT_EQ(OldMesh.Cell1DsExtrema(0, originalNumEdges), 4);
    EXPECT_EQ(OldMesh.Cell1DsExtrema(1, originalNumEdges), 5);
}

/*
TEST(TestChecking, TestCheckAddEdges2)
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

TEST(CheckAddEdgesTest, TestCheckAddEdges3) {
	
	//aggiungo un nuovo lato che non esiste
    PolyhedralMesh mesh;
    mesh.Cell1DsExtrema = MatrixXi::Zero(2,10);
	mesh.Cell1DsExtrema.col(0) = Vector2i(0,1);
	mesh.Cell1DsExtrema.col(1) = Vector2i(1,2);
	mesh.Cell1DsExtrema.col(2) = Vector2i(2,0);
	
	mesh.Cell1DsId = {0,1,2};
	
	int idE_origin = 2;
	
    Vector2i edge(2,3);

    int idE = CheckAddEdges(mesh, edge, idE_origin);

    EXPECT_EQ(idE, 3);              
}







/*TEST(TestChecking, TestCheckAddVertices1)
{
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	Vector3d vertex={0,0,1};
	int id_prova=CheckAddVertices(OldMesh,vertex,OldMesh.NumCell0Ds);
	int id_corretto=OldMesh.NumCell0Ds;
	EXPECT_EQ(id_ptova,id_corretto);
}

TEST(TestChecking, TestCheckAddVertices2)
{
	PolyhedralMesh OldMesh;
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	Vector3d vertex={0,1,1};
	int id_prova=CheckAddVertices(OldMesh,vertex,OldMesh.NumCell0Ds);
	int id_corretto=OldMesh.NumCell0Ds+1;
	EXPECT_EQ(id_ptova,id_corretto);
}

TEST(TestChecking, TestTriangulationTypeI)
{
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh NewMesh;
	int p=4;
	int q=3;
	int n=0;
	TriangulationTypeI(OldMesh,NewMesh,p,q,n);
	EXPECT_EQ();
}

TEST(TestChecking, TestTriangulationTypeII)
{
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh NewMesh;
	int n=0;
	TriangulationTypeII(OldMesh,NewMesh,n);
	EXPECT_EQ();
}

TEST(TestChecking, TestTriangulateFaces)
{
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh NewMesh;
	int p=3;
	int q=4;
	int b=2;
	int c=2;
	TriangulateFaces(OldMesh,NewMesh,p,q,b,c);
	EXPECT_EQ();
}

TEST(TestChecking, TestDualConstructor)
{
	PolyhedralMesh OldMesh;
	string file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	string file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	string file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	string file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds);
	PolyhedralMesh DualMesh;
	DualConstructor(OldMesh,DualMesh);
	EXPECT_EQ();
}*/
}