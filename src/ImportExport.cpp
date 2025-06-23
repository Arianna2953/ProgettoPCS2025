#include "ImportExport.hpp"
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
		//cout << "Cell0D ShortPath:" << endl;
		if (mesh.ShortPathCell0Ds.size()!=0){
			for(auto it = mesh.ShortPathCell0Ds.begin(); it != mesh.ShortPathCell0Ds.end(); it++){
				//cout << "ShortPath id:\t" << it -> first << "\t valori:";
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
		else{
			//cout<< "No not-null ShortPath found" << endl;
		}
	}

    if(!ImportCell1Ds(file1Ds,mesh)){
        return false;
	}
	else{
		//cout << "Cell1D ShortPath:" << endl;
		if (mesh.ShortPathCell1Ds.size() != 0){
			for(auto it = mesh.ShortPathCell1Ds.begin(); it != mesh.ShortPathCell1Ds.end(); it++){
				//cout << "ShortPath id:\t" << it -> first << "\t valori:";
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
		else{
			//cout<< "No not-null ShortPath found" << endl;
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

bool ImportCell0Ds(const string& file0Ds, PolyhedralMesh& mesh)
{
    ifstream file(file0Ds);

    if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
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
        char delimiter; // per memorizzare il ; del file csv

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
                // mesh.ShortPathCell0Ds[ShortPath].push_back(id);
                it->second.push_back(id);
            }
        }
    }
    return true;
}

bool ImportCell1Ds(const string& file1Ds, PolyhedralMesh& mesh)
{
    ifstream file(file1Ds);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
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
                // mesh.ShortPathCell1Ds[ShortPath].push_back(id);
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

bool ImportCell2Ds(const string& file1Ds, PolyhedralMesh& mesh)
{
    ifstream file(file1Ds);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
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
	//int ShortPath;
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
		
		/*if(ShortPath != 0)
        {
            const auto it = mesh.ShortPathCell2Ds.find(ShortPath);
            if(it == mesh.ShortPathCell2Ds.end())
            {
                mesh.ShortPathCell2Ds.insert({ShortPath, {id}});
            }
            else
            {
                // mesh.ShortPathCell2Ds[ShortPath].push_back(id);
                it->second.push_back(id);
            }
		}*/
		
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
			//cout << "area poligono " << id << " = " << area << endl;
		
			if (area<=1.e-16)
			{
			cerr<<"il poligono ha area pari a 0";
			return false;
			}
			
		}

    return true;
}

bool ImportCell3Ds(const string& file3Ds, PolyhedralMesh& mesh)
{
    ifstream file(file3Ds);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
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
	//int marker;
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

		
	/*if(marker != 0)
        {
            const auto it = mesh.MarkerCell3Ds.find(marker);
            if(it == mesh.MarkerCell3Ds.end())
            {
                mesh.MarkerCell3Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell3Ds[marker].push_back(id);
                it->second.push_back(id);
            }
		}*/
		
		// verifica che tutti i poligoni abbiano volume diverso da zero
		
    }
    return true;
}

bool Exportfile0Ds(const PolyhedralMesh& polyNew)
{
	string output = "Cell0Ds.txt"; //creo il file di output
	ofstream ofs(output);
	if (!ofs){ // controllo che il file di output si apra correttamente
		cerr << "Errore nell'apertura del file" << endl;
		return false;
	}
	ofs << "id,ShortPath,x,y,z" << endl;//header
	for (int i = 0; i < polyNew.NumCell0Ds; i++){
		int id = polyNew.Cell0DsId[i];		
		/*int sp;
		if (polyNew.ShortPathCell0Ds.contains(1)){
			// const l=polyNew.ShortPathCell0Ds.get(1)
			int nome=1;
			for (elemento:: polyNew.ShortPathCell0Ds.at(nome)){
				if (elemento==id){
					sp=1;
				}
				else{
					sp=0;
				}
			};
		}
		else{
			cout << "ShortPath not found" << endl;
		}*/
		//list<int> marker=polyNew.MarkerCell0Ds.at(id);
		
		int shortpathCode = -1;
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

bool Exportfile1Ds(const PolyhedralMesh & polyNew)
{
	string output="Cell1Ds.txt"; //creo il file di output
	ofstream ofs(output);
	if (!ofs){ // controllo che il file di output si apra correttamente
		cerr << "Errore nell'apertura del file" << endl;
		return false;
	}
	ofs << "id,ShortPath,v0,v1" << endl;//header
	for (int i = 0; i < polyNew.NumCell1Ds; i++){
		int id = polyNew.Cell1DsId[i];
		/*int sp;
		if (polyNew.ShortPathCell1Ds.contains(1)){
			// const l=polyNew.ShortPathCell1Ds.get(1);
			int nome=1;
			for (elemento:: polyNew.ShortPathCell1Ds.at(nome)){
				if (elemento==id){
					sp=1;
				}
				else{
					sp=0;
				}
			};
		}
			cout << "ShortPath not found" << endl;
		}*/
		//list<int> marker=polyNew.MarkerCell1Ds.at(id);
		
		int shortpathCode = -1;
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

bool Exportfile2Ds(const PolyhedralMesh & polyNew)
{
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

bool Exportfile3Ds(const PolyhedralMesh & polyNew)
{
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

bool ExportPolyhedron(const PolyhedralMesh& polyNew)
{
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
		
