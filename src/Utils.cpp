#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"
#include <string>
#include <list>
#include <map>
#include <vector>
#include <set>

namespace PolyhedralLibrary
{
bool ImportMesh(PolyhedralMesh& mesh)
{

    if(!ImportCell0Ds(mesh)){
        return false;
	}
	else{
		cout << "Cell0D marker:" << endl;
		if (mesh.MarkerCell0Ds.size()!=0){
			for(auto it=mesh.MarkerCell0Ds.begin(); it!=mesh.MarkerCell0Ds.end(); it++){
				cout << "marker id:\t" << it -> first << "\t valori:";
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
		else{
			cout<< "No not-null marker found" << endl;
		}
	}

    if(!ImportCell1Ds(mesh)){
        return false;
	}
	else{
		cout << "Cell1D marker:" << endl;
		if (mesh.MarkerCell1Ds.size()!=0){
			for(auto it=mesh.MarkerCell1Ds.begin(); it!=mesh.MarkerCell1Ds.end(); it++){
				cout << "marker id:\t" << it -> first << "\t valori:";
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
		else{
			cout<< "No not-null marker found" << endl;
		}
	}

    if(!ImportCell2Ds(mesh)){
        return false;
	}
	else{
		cout << "Cell2D marker:" << endl;
		if (mesh.MarkerCell2Ds.size()!=0){
			for(auto it=mesh.MarkerCell2Ds.begin(); it!=mesh.MarkerCell2Ds.end(); it++){
				cout << "marker id:\t" << it -> first << "\t valori:";
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
		else{
			cout<< "No not-null marker found" << endl;
		}
	}
	
	if(!ImportCell3Ds(mesh)){
        return false;
	}
	else{
		cout << "Cell3D marker:" << endl;
		if (mesh.MarkerCell3Ds.size()!=0){
			for(auto it=mesh.MarkerCell3Ds.begin(); it!=mesh.MarkerCell3Ds.end(); it++){
				cout << "marker id:\t" << it -> first << "\t valori:";
				for (const unsigned int id: it -> second)
					cout << "\t" << id;
				cout << endl;
			}
		}
		else{
			cout<< "No not-null marker found" << endl;
		}
	}

    return true;
}

bool ImportCell0Ds(& nome_file, PolygonalMesh& mesh)
{
    ifstream file("./Cell0Ds_" << nome_file << ".csv");

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

        unsigned int id;
        unsigned int marker;
        char delimiter; // per memorizzare il ; del file csv

        converter >>  id >> delimiter >> marker >> delimiter >> mesh.Cell0DsCoordinates(0, id) >> delimiter >> mesh.Cell0DsCoordinates(1, id);

        mesh.Cell0DsId.push_back(id);

        // Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell0Ds.find(marker);
            if(it == mesh.MarkerCell0Ds.end())
            {
                mesh.MarkerCell0Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell0Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }
    }
    return true;
}

bool ImportCell1Ds(& nome_file, PolygonalMesh& mesh)
{
    ifstream file("./Cell1Ds_" << nome_file << ".csv");

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

        unsigned int id;
        unsigned int marker;
        char delimiter;

        converter >>  id >> delimiter >> marker >> delimiter >>  mesh.Cell1DsExtrema(0, id) >> delimiter >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);

        // Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell1Ds.find(marker);
            if(it == mesh.MarkerCell1Ds.end())
            {
                mesh.MarkerCell1Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell1Ds[marker].push_back(id);
                it->second.push_back(id);
            }
		}
		
		// verifica che ciascun lato non abbia lunghezza zero
		int& origine=mesh.Cell1DsExtrema(0,id);
		int& fine=mesh.Cell1DsExtrema(1,id);
		if (origine==fine)
		{
			cerr<<"Il lato ha lunghezza pari a 0";
			return false;
		}
    }
    return true;
}

bool ImportCell2Ds(& nome_file, PolygonalMesh& mesh)
{
    ifstream file;
    file.open("./Cell2Ds_" << nome_file << ".csv");

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

        unsigned int id;
		unsigned int marker;
		unsigned int num_vert;
		unsigned int num_edges;
		char delimiter;
        
		converter >> id >> delimiter >> marker >> delimiter >> num_vert;

        vector<unsigned int> vecv;
		vecv.reserve(num_vert);
        for(unsigned int i = 0; i < num_vert; i++)
		{
			unsigned int vert;
            converter >> delimiter >> vert;
			vecv.push_back(vert);
		}
		mesh.Cell2DsVertices.push_back(vecv);
		
		converter >> id >> delimiter >> marker >> delimiter >> num_edges;
		
		vector<unsigned int> vece;
		vece.reserve(num_edges);
        for(unsigned int i = 0; i < num_edges; i++)
		{
			unsigned int edge;
            converter >> delimiter >> edge;
			vece.push_back(edge);
		}
		mesh.Cell2DsEdges.push_back(vece);
		mesh.Cell2DsId.push_back(id);
		
		if(marker != 0)
        {
            const auto it = mesh.MarkerCell2Ds.find(marker);
            if(it == mesh.MarkerCell2Ds.end())
            {
                mesh.MarkerCell2Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell2Ds[marker].push_back(id);
                it->second.push_back(id);
            }
		}
		// verifica che tutti i poligoni abbiano area diversa da zero
		vector<unsigned int>& vecvet=mesh.Cell2DsVertices[id];
		const unsigned int n=vecvet.size();
		double area=0.0;
		for(unsigned int i=0;i<n;i++)
		{
			const unsigned int viid=vecvet[i];
			const unsigned int vjid=vecvet[(i+1)%n]; // serve per chiudere il poligono e connettere primo e ultimo vertice
			const MatrixXd coord=mesh.Cell0DsCoordinates;
			const double Xvi=coord(0,viid);
			const double Yvi=coord(1,viid);
			const double Zvi=coord(2,viid);
			const double Xvj=coord(0,vjid);
			const double Yvj=coord(1,vjid);
			const double Zvi=coord(2,vjid);
			area=+(Xvi*Yvj)-(Xvj*Yvi);
		}
		area =abs(area/2.0);
		if (area<=1.e-16)
		{
			cerr<<"il poligono ha area pari a 0";
			return false;
		}
    }
    return true;
}
}

bool ImportCell3Ds(& nome_file, PolygonalMesh& mesh)
{
    ifstream file;
    file.open("./Cell3Ds_" << nome_file << ".csv");

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
    mesh.Cell3DsFaces.reserve(mesh.NumCell3Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
		unsigned int marker;
		unsigned int num_fac;
		char delimiter;
        
		converter >> id >> delimiter >> marker >> delimiter >> num_fac;

		mesh.Cell3DsId.push_back(id);
		
		if(marker != 0)
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
		}
		// verifica che tutti i poligoni abbiano volume diverso da zero
		
    }
    return true;
}
}
