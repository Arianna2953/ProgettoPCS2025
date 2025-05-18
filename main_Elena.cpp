
#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"
#include "ImportExport.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include <vector>

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

bool DualConstructor(const PolyhedralMesh& polyhedron, PolyhedralMesh& dual) 
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
	
	const unsigned int n = dual.NumCell1Ds*2/dual.NumCell2Ds; //numero di vertici/lati delle facce del duale
	
	cout << dual.NumCell0Ds << "  " << dual.NumCell1Ds << "  " << dual.NumCell2Ds << "  " << dual.NumCell3Ds << endl;
	  
	//Iterando sulle facce del poliedro originale mi trovo i baricentri, cioè i vertici del poliedro duale
	for (unsigned int f : polyhedron.Cell2DsId)
		{
			VectorXd barCoords(3); //Contenitore dove salverò le coordinate del baricentro
			
			//Ricavo id e coordinate dei vertici della faccia
			const vector<unsigned int> faceVertices = polyhedron.Cell2DsVertices[f];
			//o
			//const vector<unsigned int>& faceVertices = polyhedron.Cell2DsVertices[f];
			Vector3d v0;
			v0 = polyhedron.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d v1;
			v1 = polyhedron.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d v2;
			v2 = polyhedron.Cell0DsCoordinates.col(faceVertices[2]);
			
			//Calcolo le coordinate del baricentro della faccia
			barCoords = (v0+v1+v2)/3;
			
			//Proietto il punto sulla sfera di raggio 1
			double bn = barCoords.norm();
			if (bn < 1e-16) {
				cerr << "Warning: il vettore considerato ha lunghezza nulla";
				return false;
			};
			barCoords /= bn;
			
			//Aggiungo il punto appena trovato ai vertici del  poliedro duale
			dual.Cell0DsId.push_back(f);
			dual.Cell0DsCoordinates.col(f) << barCoords;
		};
	
	
	
	unsigned int newEdge = 0;
	for (unsigned int f : polyhedron.Cell2DsId)
	{
		for (unsigned int e : polyhedron.Cell2DsEdges[f]) 
		{
			for (unsigned int g : polyhedron.Cell2DsId)
			{
				for (unsigned int a : polyhedron.Cell2DsEdges[g])
				{
					if (e == a and f > g) 
					{
						dual.Cell1DsId.push_back(newEdge);
						dual.Cell1DsExtrema(0,newEdge) = f;
						dual.Cell1DsExtrema(1,newEdge) = g;
						

						
						newEdge++;
					}
				}
			}
		}
	}
	
	cout << dual.Cell1DsExtrema << endl;
	
	//Iterando sui vertici del  poliedro originale costruisco le facce del duale
	for (unsigned int v : polyhedron.Cell0DsId)
	{
		vector<unsigned int> adjacentFaces;
		adjacentFaces.reserve(n);
		unsigned int i = 0;
		//Trovo le vacce adiacenti al vertice v (nel poliedro originale)
		for  (unsigned int f : polyhedron.Cell2DsId)
		{
			for (unsigned int vf : polyhedron.Cell2DsVertices[f])
			{
				if (v == vf) 
				{
					adjacentFaces.push_back(f);
					i++;
					break;
				}
			}
			if (i >= n) {break;};
		}
		dual.Cell2DsId.push_back(v);
		dual.Cell2DsVertices.push_back(adjacentFaces);
		
		
		vector<unsigned int> faceEdges;
		faceEdges.reserve(n);
		for (auto it1 = adjacentFaces.begin(); it1 != adjacentFaces.end(); it1++)
		{
			for (auto it2 = adjacentFaces.begin(); it2 < it1 and it2 != adjacentFaces.end(); it2++)
			{
				Vector2i ex1;
				ex1 << *it1, *it2;
				Vector2i ex2;;
				ex2 << *it2, *it1;
				for (unsigned int j = 0; j < dual.NumCell1Ds; j++)
				{
					if (ex1 == dual.Cell1DsExtrema.col(j) or ex2 == dual.Cell1DsExtrema.col(j))
					{
						faceEdges.push_back(j);
					};
				}
			}
		}
		dual.Cell2DsEdges.push_back(faceEdges);
		
		
	}
	
return true;	
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
	
	/*int p = 3;
    int q = 3;
    int b = 2;
    int c = 0;
    int v0 = 0, v1 = 0;*/
	
	PolyhedralMesh poly1;
	PolyhedralMesh poly2;
	
	if(!ImportMesh(poly1,file0Ds,file1Ds,file2Ds,file3Ds))
	{
		cerr << "File non trovato." << endl;
		return 1;
	};
	
	
	if(!DualConstructor(poly1, poly2))
	{
		cerr << "Costruzione del duale non riuscita." << endl;
		return 1;
	};
	
	return 0;
}

/*
1. calcola baricentro di una faccia
		sia n numero di facce
		matrice nx5 in cui memorizzo id_faccia, id_baricentro, coordinate_baricentro
	
		
2. connetti baricentri

3. 
*/
