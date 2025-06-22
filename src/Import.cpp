#include "Import.hpp"
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
bool ImportMesh(PolyhedralMesh& mesh,const string& file0Ds,const string& file1Ds,const string& file2Ds,const string& file3Ds){

    if(!ImportCell0Ds(file0Ds,mesh)){
        return false;
	}
	else{
		if (mesh.ShortPathCell0Ds.size()!=0){
			for(auto it = mesh.ShortPathCell0Ds.begin(); it != mesh.ShortPathCell0Ds.end(); it++){
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
	}

    if(!ImportCell1Ds(file1Ds,mesh)){
        return false;
	}
	else{
		if (mesh.ShortPathCell1Ds.size() != 0){
			for(auto it = mesh.ShortPathCell1Ds.begin(); it != mesh.ShortPathCell1Ds.end(); it++){
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
	}

    if(!ImportCell2Ds(file2Ds,mesh)){
        return false;
	}
	
	if(!ImportCell3Ds(file3Ds,mesh)){
        return false;
	}

    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool ImportCell0Ds(const string& file0Ds, PolyhedralMesh& mesh){
    
	ifstream file(file0Ds);
	if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0) //controllo che il file non sia vuoto
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        int id;
        int ShortPath;
        char delimiter; // per memorizzare il ";" del file csv

        converter >>  id >> delimiter >> ShortPath >> delimiter >> mesh.Cell0DsCoordinates(0, id) >> delimiter >> mesh.Cell0DsCoordinates(1, id)>> delimiter >> mesh.Cell0DsCoordinates(2, id);

        mesh.Cell0DsId.push_back(id);

        // Memorizza gli ShortPath
        if(ShortPath != 0)
        {
            const auto it = mesh.ShortPathCell0Ds.find(ShortPath);
            if(it == mesh.ShortPathCell0Ds.end())
            {
                mesh.ShortPathCell0Ds.insert({ShortPath, {id}});
            }
            else
            {
                it->second.push_back(id);
            }
        }
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool ImportCell1Ds(const string& file1Ds, PolyhedralMesh& mesh){
    
	ifstream file(file1Ds);
    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, mesh.NumCell1Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        int id;
        int ShortPath;
        char delimiter;

        converter >>  id >> delimiter >> ShortPath >> delimiter >>  mesh.Cell1DsExtrema(0, id) >> delimiter >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);

        // Memorizza gli ShortPath
        if(ShortPath != 0)
        {
            const auto it = mesh.ShortPathCell1Ds.find(ShortPath);
            if(it == mesh.ShortPathCell1Ds.end())
            {
                mesh.ShortPathCell1Ds.insert({ShortPath, {id}});
            }
            else
            {
                it->second.push_back(id);
            }
		}
		
		// verifica che ciascun lato non abbia lunghezza zero
		if(mesh.Cell1DsExtrema(0, id) == mesh.Cell1DsExtrema(1, id)){
			cerr<<"Il lato ha lunghezza pari a 0";
			return false;
		}
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool ImportCell2Ds(const string& file1Ds, PolyhedralMesh& mesh){
    
	ifstream file(file1Ds);
    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo header
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);
		int id;
		int num_vert;
		int num_edges;
		char delimiter;
		
		converter >> id >> delimiter >> num_vert;

        vector<int> vecv;
		vecv.reserve(num_vert);
        for(int i = 0; i < num_vert; i++){
			int vert;
			converter >> delimiter >> vert;
			vecv.push_back(vert);
		}
		mesh.Cell2DsVertices.push_back(vecv);
		
		converter >> delimiter >> num_edges;
		
		vector<int> vece;
		vece.reserve(num_edges);
        for(int i = 0; i < num_edges; i++)
		{
			int edge;
            converter >> delimiter >> edge;
			vece.push_back(edge);
		}
		mesh.Cell2DsEdges.push_back(vece);
		mesh.Cell2DsId.push_back(id);
		
		// verifica che tutti i poligoni abbiano area diversa da zero
		const MatrixXd& coord = mesh.Cell0DsCoordinates;
			Vector3d areaVec(0.0, 0.0, 0.0);

			for (int i = 0; i < num_vert; i++) {
				int j = (i + 1) % num_vert;

				Vector3d vi = coord.col(vecv[i]);
				Vector3d vj = coord.col(vecv[j]);
				areaVec += vi.cross(vj);
			}

			double area = 0.5 * areaVec.norm();
			if (area<=1.e-16)
			{
			cerr<<"il poligono ha area pari a 0";
			return false;
			}
			
		}

    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool ImportCell3Ds(const string& file3Ds, PolyhedralMesh& mesh){
    
	ifstream file(file3Ds);
	if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    //rimuovo header
    listLines.pop_front();

    mesh.NumCell3Ds = listLines.size();

    if (mesh.NumCell3Ds == 0)
    {
        cerr << "There is no cell 3D" << endl;
        return false;
    }
	
	mesh.Cell3DsId.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsVertices.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsEdges.reserve(mesh.NumCell3Ds);
	mesh.Cell3DsFaces.reserve(mesh.NumCell3Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        int id;
		int num_vert;
		int num_edges;
		int num_faces;
		char delimiter;
        
		converter >> id >> delimiter >> num_vert;

        vector<int> vecv;
		vecv.reserve(num_vert);
        for(int i = 0; i < num_vert; i++){
			int vert;
			converter >> delimiter >> vert;
			vecv.push_back(vert);
		}
		mesh.Cell3DsVertices.push_back(vecv);
		
		converter >> delimiter >> num_edges;
		
		vector<int> vece;
		vece.reserve(num_edges);
        for(int i = 0; i < num_edges; i++)
		{
			int edge;
            converter >> delimiter >> edge;
			vece.push_back(edge);
		}
		mesh.Cell3DsEdges.push_back(vece);
		mesh.Cell3DsId.push_back(id);
     
		converter >> delimiter >> num_faces;
	
		vector<int> vecf;
		vecf.reserve(num_faces);
		for(int i = 0; i < num_faces; i++)
		{
			int face;
            converter >> delimiter >> face;
			vecf.push_back(face);
		}
		mesh.Cell3DsFaces.push_back(vecf);		
    }
    return true;
}

}

