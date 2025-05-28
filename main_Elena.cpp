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

bool Du1alConstructor(const PolyhedralMesh& polyhedron, PolyhedralMesh& dual) 
{	
	const int V = polyhedron.NumCell2Ds; //Numero facce poliedro originale = numero vertici duale
	const int E = polyhedron.NumCell1Ds; //Numero lati poliedro originale = numero lati duale
	const int F = polyhedron.NumCell0Ds; //Numero vertici poliedro originale = numero facce duale
	
	dual.NumCell0Ds = V;
	dual.Cell0DsId.reserve(V); //Vettore id vertici del duale
	dual.Cell0DsCoordinates = MatrixXd::Zero(3,V); //Matrice coordinate vertici del duale
	
	dual.NumCell1Ds = E;
	dual.Cell1DsId.reserve(E); //Vettore id lati del duale
	dual.Cell1DsExtrema = MatrixXi::Zero(2,E); //Matrice id estremi dei lati del duale
	
	dual.NumCell2Ds = F;
	dual.Cell2DsId.reserve(F); //Vettore id facce del duale
	dual.Cell2DsVertices.reserve(F);
	dual.Cell2DsEdges.reserve(F);
	
	dual.NumCell3Ds = 1;
	dual.Cell3DsId.push_back(0);
	dual.Cell3DsEdges.reserve(1);
	dual.Cell3DsVertices.reserve(1);
	dual.Cell3DsFaces.reserve(1);
	
	const int Npoly = 3; //Numero di vertici/lati per faccia nel poliedro originale 
	const int Ndual = E*2/F; //Numero di vertici/lati per faccia nel poliedro duale 
	
	//Iterando sulle facce del poliedro originale mi trovo i baricentri, cioè i vertici del poliedro duale
	for (int f : polyhedron.Cell2DsId)
		{
			Vector3d barCoords; //Contenitore dove salverò le coordinate del baricentro
			
			//Ricavo id e coordinate dei vertici della faccia
			const vector<int>& faceVertices = polyhedron.Cell2DsVertices[f];
			Vector3d v0 = polyhedron.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d v1 = polyhedron.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d v2 = polyhedron.Cell0DsCoordinates.col(faceVertices[2]);
			
			//Calcolo le coordinate del baricentro della faccia
			barCoords = (v0+v1+v2)/3;
			
			//Proietto il punto sulla sfera di raggio 1
			if (barCoords.norm() < 1e-16) {
				cerr << "Warning: il vettore considerato ha lunghezza nulla";
				return false;
			};
			barCoords.normalize();
			
			//Aggiungo il punto appena trovato ai vertici del  poliedro duale
			dual.Cell0DsId.push_back(f);
			dual.Cell0DsCoordinates.col(f) << barCoords;
		};
	
	/*Trovo i lati del duale. 
	Scorro le facce (f1) del poliedro orginale e ne considero un lato alla volta cercando la faccia (f2) 
	con cui è condiviso. I baricentri di f1 e f2, nel duale, saranno gli estremi di un lato che vado ad aggiungere all'elenco. */
	int newEdge = 0; //Id lato da aggiungere
	for (int i = 0; i < polyhedron.NumCell2Ds; i++) //Itero sulle facce del poliedro originale
	{
		int f1 = polyhedron.Cell2DsId[i];
		vector<int> edgesf1 = polyhedron.Cell2DsEdges[f1];
		for (int h = 0; h < Npoly; h++) //Considero 1 alla volta i lati della faccia f1 
		{
			bool found = false; //Variabile booleana indica se ho trovato o no "l'altra faccia" del lato edgesf1[h]. NOTA: ogni lato è condiviso solo da 2 facce.
			for (int j = 0; j <i; j ++) //Considero una alla volta le facce del poliedro originale "precedenti" a f1
			{
				int f2 = polyhedron.Cell2DsId[j];
				vector<int> edgesf2 = polyhedron.Cell2DsEdges[f2];
				//bool found = false; //Variabile booleana indica se ho trovato o no "l'altra faccia" del lato edgesf1[h]. NOTA: ogni lato è condiviso solo da 2 facce.
				for (int k = 0; k < Npoly; k++)
				{
					if (edgesf1[h] == edgesf2[k]){
						dual.Cell1DsId.push_back(newEdge);
						dual.Cell1DsExtrema(0,newEdge) = f1;
						dual.Cell1DsExtrema(1,newEdge) = f2;
						newEdge++; //Incremento id nuovo lato.
						found = true; //Ho trovato faccia in comune.
						break; //Ho trovato edgesf1[h] tra i lati di f2 non ha senso considerare gli altri lati di f2.
					} 
				}
				if (found) {break;} //Ho trovato la "seconda faccia" di edgesf1[h], posso passare a condiderare edgesf1[h+1].
			}
		}
	}
	
	
	//Iterando sui vertici del  poliedro originale costruisco le facce del duale
	for (int v : polyhedron.Cell0DsId)
	{
		vector<int> adjacentFaces; //Contenitore dove mi salverò gli id delle facce adiacenti a v nel poliedro originale.
		adjacentFaces.reserve(Ndual);
		int m = 0; //Contatore delle facce adiacenti a v trovate
		//Trovo le vacce adiacenti al vertice v (nel poliedro originale)
		for  (int f : polyhedron.Cell2DsId)
		{
			for (int vf : polyhedron.Cell2DsVertices[f])
			{
				if (v == vf) 
				{
					adjacentFaces.push_back(f);
					m++;
					break; //OSS: il vertice v non può appartenere alla faccia f più di una volta. se l'ho trovato tra i suoi vertici posso passare alla faccia successiva.
				}
			}
			if (m >= Ndual) {break;}; //Ci sono al più Ndual facce adiacenti a v.
		}
		vector<int>& faceVerticesDual = adjacentFaces; //Quelli che nel poliedro originale erano gli id delle facce adiacenti al vertice v, nel duale sono gli id dei vertici della faccia con id == v.
		
		vector<int> faceEdgesDual; 
		faceEdgesDual.reserve(Ndual);
		vector<int> fVDSorted;
		fVDSorted.reserve(Ndual);
		
		auto it0 = faceVerticesDual.begin();
		auto it1 = faceVerticesDual.begin();
		for (int NoE = 0; NoE < Ndual; NoE++)
		{
			bool foundE = false;
			for (auto it2 = faceVerticesDual.begin(); it2 != faceVerticesDual.end(); it2++)
			{	
				if (it2!=it1 && it2!=it0) 
				{
					Vector2i ex1;
					ex1 << *it1, *it2;
					Vector2i ex2;;
					ex2 << *it2, *it1;
					for (int j = 0; j < E; j++) //Scorro i lati del duale
					{
						if (ex1 == dual.Cell1DsExtrema.col(j) || ex2 == dual.Cell1DsExtrema.col(j))
						{
							faceEdgesDual.push_back(dual.Cell1DsId[j]);
							fVDSorted.push_back(*it1);
							it0 = it1;
							it1 = it2;
							foundE = true;
							break;
						}
					}
				}
				if (foundE) {break;}	
			}
		}
		
		dual.Cell2DsId.push_back(v); //v corrisponde ad una faccia del duale.
		dual.Cell2DsVertices.push_back(fVDSorted); 
		dual.Cell2DsEdges.push_back(faceEdgesDual);
		
	}
	
	
	//Aggiorno valori Cell3Ds
	dual.Cell3DsVertices.push_back(dual.Cell0DsId);
	dual.Cell3DsEdges.push_back(dual.Cell1DsId);
	dual.Cell3DsFaces.push_back(dual.Cell2DsId);
	return true;	
};


int main() 
{
	string file0Ds;
    string file1Ds;
    string file2Ds;
    string file3Ds;
	
	file0Ds = "./Cell0Ds_octahedron.csv";
	file1Ds = "./Cell1Ds_octahedron.csv";
	file2Ds = "./Cell2Ds_octahedron.csv";
	file3Ds = "./Cell3Ds_octahedron.csv";
	
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
	
	if (!ExportPolyhedron(poly2))
	{
		cerr << "Esportazione fallita." << endl;
		return 1;
	}
	
	Gedim::UCDUtilities utilities;
    {	vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(poly2.NumCell0Ds, 0.0);
        for(const auto &m : poly2.MarkerCell0Ds)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();

        utilities.ExportPoints("./Cell0Ds.inp",
                               poly2.Cell0DsCoordinates,
                               cell0Ds_properties);
    }

    {	vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(poly2.NumCell1Ds, 0.0);
        for(const auto &m : poly2.MarkerCell1Ds)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 poly2.Cell0DsCoordinates,
                                 poly2.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
	}
	
	return 0;
}