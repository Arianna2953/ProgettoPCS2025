#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"


using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

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
	
	PolyhedralMesh poly1;
	
	if(!ImportMesh(poly1,file0Ds,file1Ds,file2Ds,file3Ds))
		{
			cerr << "file non trovato" << endl;
			return 1;
		}
	
	
	return 0;
}

/*
1. calcola baricentro di una faccia
		sia n numero di facce
		matrice nx5 in cui memorizzo id_faccia, id_baricentro, coordinate_baricentro
		
2. connetti baricentri

3. 
*/