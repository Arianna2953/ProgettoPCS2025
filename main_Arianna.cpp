#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include <cmath> 
#include <numeric> 
#include <chrono>
#include <vector>
#include "Import.hpp"
#include "Export.hpp"
#include "Triangulation.hpp"
#include "Dual.hpp"
#include "ShortPath.hpp"
#include "PolyhedralMesh.hpp"
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
	
	file0Ds = "./Cell0Ds_icosahedron.csv";
	file1Ds = "./Cell1Ds_icosahedron.csv";
	file2Ds = "./Cell2Ds_icosahedron.csv";
	file3Ds = "./Cell3Ds_icosahedron.csv";
	
	int p = 3;
    int q = 3;
    int b = 2;
    int c = 0;
    int v0 = 0, v1 = 0;
	
	PolyhedralMesh poly1;
	PolyhedralMesh poly2;
	
	if(!ImportMesh(poly1,file0Ds,file1Ds,file2Ds,file3Ds))
	{
		cerr << "File non trovato." << endl;
		return 1;
	};
	if(!ExportPolyhedron(poly1))
    {
        cerr << "file not found" << endl;
        return false;
    }
	
}