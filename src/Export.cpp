#include "Export.hpp"
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

using namespace std;

namespace PolyhedralLibrary
{
	
bool Exportfile0Ds(const PolyhedralMesh& polyNew){
	string output = "Cell0Ds.txt"; //creo il file di output
	ofstream ofs(output);
	if (!ofs){ // controllo che il file di output si apra correttamente
		cerr << "Errore nell'apertura del file" << endl;
		return false;
	}
	ofs << "id,ShortPath,x,y,z" << endl;//header
	for (int i = 0; i < polyNew.NumCell0Ds; i++){
		int id = polyNew.Cell0DsId[i];				
		int shortpathCode = 0;
        for (const auto& [code, ids] : polyNew.ShortPathCell0Ds) {
            if (find(ids.begin(), ids.end(), id) != ids.end()) {
                shortpathCode = code;
                break;
            }
        }
        ofs << id << "," << shortpathCode << ",";
		ofs << polyNew.Cell0DsCoordinates(0, id) << ";" << polyNew.Cell0DsCoordinates(1, id) << ";" << polyNew.Cell0DsCoordinates(2, id) << endl;
	}
		
	ofs.close();
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool Exportfile1Ds(const PolyhedralMesh & polyNew){
	string output="Cell1Ds.txt"; //creo il file di output
	ofstream ofs(output);
	if (!ofs){ // controllo che il file di output si apra correttamente
		cerr << "Errore nell'apertura del file" << endl;
		return false;
	}
	ofs << "id,ShortPath,v0,v1" << endl;//header
	for (int i = 0; i < polyNew.NumCell1Ds; i++){
		int id = polyNew.Cell1DsId[i];
		
		int shortpathCode = 0;
		for (const auto& [code, ids] : polyNew.ShortPathCell1Ds) {
			if (find(ids.begin(), ids.end(), id) != ids.end()) {
				shortpathCode = code;
				break;
			}
		}
		ofs << id << "," << shortpathCode << ",";
		ofs << polyNew.Cell1DsExtrema(0, id) << ";" << polyNew.Cell1DsExtrema(1, id) << endl;
	}	
	
	ofs.close();
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool Exportfile2Ds(const PolyhedralMesh & polyNew){
	string output = "Cell2Ds.txt"; //creo il file di output
	ofstream ofs(output);
	if (!ofs){ // controllo che il file di output si apra correttamente
		cerr << "Errore nell'apertura del file" << endl;
		return false;
	}
	ofs << "id,n_vertices,vertices,n_edges,edges"<<endl;
	for (int i = 0; i<polyNew.NumCell2Ds; i++){
		int id = polyNew.Cell2DsId[i];
		ofs << id << ";" << polyNew.Cell2DsVertices[i].size() << ";";
		for (size_t k = 0; k<polyNew.Cell2DsVertices[i].size(); k++){
			ofs << polyNew.Cell2DsVertices[i][k] << ";";
		}
		ofs << polyNew.Cell2DsEdges[i].size() << ";";
		for (size_t j = 0; j < polyNew.Cell2DsEdges[i].size(); j++){
			ofs << polyNew.Cell2DsEdges[i][j] << ";";
		}
		ofs << endl;	
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool Exportfile3Ds(const PolyhedralMesh & polyNew){
	string output = "Cell3Ds.txt"; //creo il file di output
	ofstream ofs(output);
	if (!ofs){ // controllo che il file di output si apra correttamente
		cerr << "Errore nell'apertura del file" << endl;
		return false;
	}
	ofs << "id,n_vertices,vertices,n_edges,edges,n_faces,faces" << endl;
	for (int i = 0; i < polyNew.NumCell3Ds; i++){
		int id = polyNew.Cell3DsId[i];
		ofs << id << ";" << polyNew.Cell3DsVertices[i].size() << ";";
		for (size_t k = 0; k < polyNew.Cell3DsVertices[i].size(); k++){
			ofs << polyNew.Cell3DsVertices[i][k] << ";";
		}
		ofs << polyNew.Cell3DsEdges[i].size() << ";";
		for (size_t j = 0; j < polyNew.Cell3DsEdges[i].size(); j++){
			ofs << polyNew.Cell3DsEdges[i][j] << ";";
		}
		ofs << polyNew.Cell3DsFaces[i].size() << ";";
		for (size_t j = 0; j < polyNew.Cell3DsFaces[i].size(); j++){
			ofs << polyNew.Cell3DsFaces[i][j] << ";";
		}
		ofs << endl;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool ExportPolyhedron(const PolyhedralMesh& polyNew){
	if(!Exportfile0Ds(polyNew))
    {
        cerr << "file not found" << endl;
        return false;
    }
	
	if(!Exportfile1Ds(polyNew))
    {
        cerr << "file not found" << endl;
        return false;
    }
	
	if(!Exportfile2Ds(polyNew))
    {
        cerr << "file not found" << endl;
        return false;
    }
	
	if(!Exportfile3Ds(polyNew))
    {
        cerr << "file not found" << endl;
        return false;
    }
	return true;
}	

}	
