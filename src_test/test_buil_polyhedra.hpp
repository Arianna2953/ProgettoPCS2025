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
	OldMesh.NumCell0Ds=6;
	OldMesh.Cell0DsId={0,1,2,3,4,5};
	OldMesh.Cell0DsCoordinates[3][mesh.NumCell0Ds]={
		{1,0,0},
		{-1,0,0},
		{0,1,0},
		{0,-1,0},
		{0,0,1},
		{0,0,-1}};
	OldMesh.ShortPathCell0Ds={0,{0,1,2,3,4,5}};
	OldMesh.NumCell1Ds=12;
	OldMesh.Cell1DsId={0,1,2,3,4,5,6,7,8,9,10,11};
	OldMesh.Cell1DsExtrema[2][mesh.NumCell1Ds]={
		{0,2},
		{0,3},
		{0,4},
		{0,5},
		{1,2},
		{1,3},
		{1,4},
		{1,5},
		{2,4},
		{2,5},
		{3,4},
		{3,5}};
	OldMesh.ShortPathCell1Ds={0,{0,1,2,3,4,5,6,7,8,9,10,11}};
	OldMesh.NumCell2Ds=8;
	OldMesh.Cell2DsId={0,1,2,3,4,5,6,7};
	OldMesh.Cell2DsVertices={{0,2,4},{2,1,4},{1,3,4},{3,0,4},{0,2,5},{2,1,5},{1,3,5},{3,0,5}};
	OldMesh.Cell2DsEdges={{0,8,2},{4,6,8},{5,10,6},{1,2,10},{0,9,3},{4,7,9},{5,11,7},{1,3,11}};
	OldMesh.NumCell3Ds=1;
	OldMesh.Cell3DsId={0};
	OldMesh.Cell3DsVertices={0,1,2,3,4,5};
	OldMesh.Cell3DsEdges={0,1,2,3,4,5,6,7,8,9,10,11};
	OldMesh.Cell3DsFaces={0,1,2,3,4,5,6,7};
	
	PolyhedralMesh NewMesh;*/
TEST(TestChecking, TestCheckAddEdges1)
{
	PolyhedralMesh OldMesh;
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
	Vector2i edge={0,5};
	int id_prova=CheckAddEdges(OldMesh,edge,OldMesh.NumCell1Ds);
	int id_corretto=OldMesh.NumCell1Ds;
	EXPECT_EQ(id_prova,id_corretto);
}

TEST(TestChecking, TestCheckAddEdges2)
{
	PolyhedralMesh OldMesh;
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
	Vector2i edge={1,5};
	int id_prova=CheckAddEdges(OldMesh,edge,OldMesh.NumCell1Ds);
	int id_corretto=OldMesh.NumCell1Ds+1;
	EXPECT_EQ(id_prova,id_corretto);
}

TEST(TestChecking, TestCheckAddVertices1)
{
	PolyhedralMesh OldMesh;
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
	Vector3d vertex={0,0,1};
	int id_prova=CheckAddVertices(OldMesh,vertex,OldMesh.NumCell0Ds);
	int id_corretto=OldMesh.NumCell0Ds;
	EXPECT_EQ(id_ptova,id_corretto);
}

TEST(TestChecking, TestCheckAddVertices2)
{
	PolyhedralMesh OldMesh;
	PolyhedralMesh OldMesh;
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
	Vector3d vertex={0,1,1};
	int id_prova=CheckAddVertices(OldMesh,vertex,OldMesh.NumCell0Ds);
	int id_corretto=OldMesh.NumCell0Ds+1;
	EXPECT_EQ(id_ptova,id_corretto);
}

TEST(TestChecking, TestTriangulationTypeI)
{
	PolyhedralMesh OldMesh;
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
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
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
	PolyhedralMesh NewMesh;
	int n=0;
	TriangulationTypeII(OldMesh,NewMesh,n);
	EXPECT_EQ();
}

TEST(TestChecking, TestTriangulateFaces)
{
	PolyhedralMesh OldMesh;
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
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
	file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
	file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
	file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
	file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
	ImportMesh(OldMesh,file0Ds,file1Ds,file2Ds,file3Ds)
	PolyhedralMesh DualMesh;
	DualConstructor(OldMesh,DualMesh);
	EXPECT_EQ();
}
}