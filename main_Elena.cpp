
#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"
#include "ImportExport.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

void DualConstructor(const PolyhedralMesh& polyhedron, PolyhedralMesh& dual) 
{	
	/*//Dati poliedro originale
	unsigned int numVPoly = polyhedron.NumCell0Ds; //numero vertici poliedro originale == numero facce del duale
	unsigned int numFPoly = polyhedron.NumCell2Ds; //numero facce poliedro originale == numero vertici del duale
	vector<unsigned int> verticesPoly = polyhedron.Cell0DsId; //id vertici poliedro originale
	vector<unsigned int> facesPoly = polyhedron.Cell1DsId; //id facce poliedro originale
	
	//Numero verici e numero facce duale
	dual.NumCell0Ds = numVPoly;
	dual.NumCell2Ds = numFPoly;
//COME CALCOLARE NUMERO LATI?	
	*/
	
	
	//Numero vertici, lati e facce del duale.
	dual.NumCell0Ds = polyhedron.NumCell2Ds;
	dual.NumCell1Ds = polyhedron.NumCell1Ds;
	dual.NumCell2Ds = polyhedron.NumCell0Ds;
	dual.NumCell3Ds = 1;
	
	dual.Cell0DsId.reserve(dual.NumCell0Ds); //Vettore id vertici del duale
	dual.Cell0DsCoordinates = MatrixXd::Zero(3,dual.NumCell0Ds); //Matrice coordinate vertici del duale
	
	dual.Cell1DsId.reserve(dual.NumCell1Ds); //Vettore id lati del duale
	dual.Cell1DsExtrema = MatrixXi::Zero(2,dual.NumCell1Ds); //Matrice id estremi dei lati del duale
	
	dual.Cell2DsId.reserve(dual.NumCell2Ds); //Vettore id facce del duale
	dual.Cell2DsVertices.reserve(dual.NumCell2Ds);
	dual.Cell2DsEdges.reserve(dual.NumCell2Ds);
	
	cout << dual.NumCell0Ds << "  " << dual.NumCell1Ds << "  " << dual.NumCell2Ds << "  " << dual.NumCell3Ds << endl;
	/*  
	//Iterando sulle facce del poliedro originale mi trovo i baricentri, cioè i vertici del poliedro duale
	for (unsigned int f in polyhedron.NumCell2Ds)
		{
			VectorXd barCoords(3); //Contenitore dove salverò le coordinate del baricentro
			
			//Ricavo id e coordinate dei vertici della faccia
			vector<unsigned int> faceVertices = polyhedron.Cell2DsVertices[f];
			VectorXd v0(3);
			v0 << polyhedron.Cell0DsCoordinates(:,faceVertices[0]);
			VectorXd v1(3);
			v1 << polyhedron.Cell0DsCoordinates(:,faceVertices[1]);
			VectorXd v2(3);
			v2 << polyhedron.Cell0DsCoordinates(:,faceVertices[2]);
			
			//Calcolo le coordinate del paricentro della faccia
			barCoords = (v0+v1+v2)/3;
		};
	*/
return;	
};


int main() 
{
	string file0Ds;
    string file1Ds;
    string file2Ds;
    string file3Ds;
	
	file0Ds = "./Cell0Ds_tetrahedron.csv";
	file1Ds = "./Cell1Ds_tetrahedron.csv";
	file2Ds = "./Cell2Ds_tetrahedron.csv";
	file3Ds = "./Cell3Ds_tetrahedron.csv";
	
	int p = 3;
    int q = 3;
    int b = 2;
    int c = 0;
    int v0 = 0, v1 = 0;
	
	PolyhedralMesh poly1;
	PolyhedralMesh poly2;
	
	if(!ImportMesh(poly1,file0Ds,file1Ds,file2Ds,file3Ds))
		{
			cerr << "file non trovato" << endl;
			return 1;
		}
	
	
	DualConstructor(poly1, poly2);
	
	return 0;
}

/*
1. calcola baricentro di una faccia
		sia n numero di facce
		matrice nx5 in cui memorizzo id_faccia, id_baricentro, coordinate_baricentro
	
		
2. connetti baricentri

3. 
