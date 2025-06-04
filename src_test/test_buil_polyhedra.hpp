#pragma once

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

namespace PolyhedralLibrary{
	PolyhedralMesh Poly;
	/*Poly.NumCell0Ds=6;
	Poly.Cell0DsId={0,1,2,3,4,5};
	Poly.Cell0DsCoordinates[3][mesh.NumCell0Ds]={
		{1,0,0},
		{-1,0,0},
		{0,1,0},
		{0,-1,0},
		{0,0,1},
		{0,0,-1}};
	Poly.ShortPathCell0Ds={0,{0,1,2,3,4,5}};
	Poly.NumCell1Ds=12;
	Poly.Cell1DsId={0,1,2,3,4,5,6,7,8,9,10,11};
	Poly.Cell1DsExtrema[2][mesh.NumCell1Ds]={
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
	Poly.ShortPathCell1Ds={0,{0,1,2,3,4,5,6,7,8,9,10,11}};
	Poly.NumCell2Ds=8;
	Poly.Cell2DsId={0,1,2,3,4,5,6,7};
	Poly.Cell2DsVertices={{0,2,4},{2,1,4},{1,3,4},{3,0,4},{0,2,5},{2,1,5},{1,3,5},{3,0,5}};
	Poly.Cell2DsEdges={{0,8,2},{4,6,8},{5,10,6},{1,2,10},{0,9,3},{4,7,9},{5,11,7},{1,3,11}};
	Poly.NumCell3Ds=1;
	Poly.Cell3DsId={0};
	Poly.Cell3DsVertices={0,1,2,3,4,5};
	Poly.Cell3DsEdges={0,1,2,3,4,5,6,7,8,9,10,11};
	Poly.Cell3DsFaces={0,1,2,3,4,5,6,7};*/
	
	PolyhedralMesh NewMesh;

TEST(TestChecking, TestCheckAddEdges)
{
	PolyhedralMesh OldMesh;
	OldMesh.Cell1DsExtrema[2][OldMesh.NumCell1Ds]={
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
	int id=5;
	CheckAddEdges(OldMesh,OldMesh.Cell1DsExtrema,id);

	EXPECT_EQ();
}

TEST(TestChecking, TestCheckAddVertices)
{
	PolyhedralMesh OldMesh;
	OldMesh.Cell0DsCoordinates[3][OldMesh.NumCell0Ds]={
		{1,0,0},
		{-1,0,0},
		{0,1,0},
		{0,-1,0},
		{0,0,1},
		{0,0,-1}};
	int id=5;
	CheckAddVertices(OldMesh,OldMesh.Cell0DsCoordinates,id);
	EXPECT_EQ();
}

TEST(TestChecking, TestTriangulationTypeI)
{
	PolyhedralMesh OldMesh;
	PolyhedralMesh NewMesh;
	TriangulationTypeI(OldMesh,NewMesh,3,3,2);
	EXPECT_EQ();
}

TEST(TestChecking, TestTriangulationTypeII)
{
	PolyhedralMesh OldMesh;
	PolyhedralMesh NewMesh;
	TriangulationTypeII(OldMesh,NewMesh,2);
	EXPECT_EQ();
}

TEST(TestChecking, TestTriangulateFaces)
{
	PolyhedralMesh OldMesh;
	TriangulateFaces(OldMesh,NewMesh,3,3,2,2);
	EXPECT_EQ();
}

TEST(TestChecking, TestDualConstructor)
{
	PolyhedralMesh OldMesh;
	PolyhedralMesh DualMesh;
	DualConstructor(OldMesh,DualMesh);
	EXPECT_EQ();
}
}