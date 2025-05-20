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
	
	const int M = E*2/V; //Numero di vertici/lati per faccia nel poliedro originale
	const int N = E*2/F; //Numero di vertici/lati per faccia nel poliedro duale
	cout << "M : " << M << "\nN : " << N << endl;
	
	//Iterando sulle facce del poliedro originale mi trovo i baricentri, cioè i vertici del poliedro duale
	for (int f : polyhedron.Cell2DsId)
		{
			VectorXd barCoords(3); //Contenitore dove salverò le coordinate del baricentro
			
			//Ricavo id e coordinate dei vertici della faccia
			const vector<int> faceVertices = polyhedron.Cell2DsVertices[f];
			//o
			//const vector<unsigned int>& faceVertices = polyhedron.Cell2DsVertices[f];
			Vector3d v0 = polyhedron.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d v1 = polyhedron.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d v2 = polyhedron.Cell0DsCoordinates.col(faceVertices[2]);
			
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
	
	cout << "Cell0DsId: ";
	for (int v : dual.Cell0DsId) {cout << v << " ";};
	cout << endl;
	cout << "Cell0DsCoordinates \n" << dual.Cell0DsCoordinates << endl;

	/*Trovo i lati del duale. 
	Scorro le facce (f1) del poliedro orginale e ne considero un lato alla volta cercando la faccia (f2) 
	con cui è condiviso. Nel duale f1 e f2 saranno gli estremi di un lato che vado ad aggiungere all'elenco. */
	
	
	int newEdge = 0; //Id lato da aggiungere
	for (int i = 0; i < polyhedron.NumCell2Ds; i++) //Itero sulle facce del poliedro originale
	{
		int f1 = polyhedron.Cell2DsId[i];
		vector<int> edges1 = polyhedron.Cell2DsEdges[f1];
		cout << edges1.size() << endl;
		for (int h = 0; h < M; h++) //Considero 1 alla volta i lati della faccia f1 
		{
			for (int j = 0; j <i; j ++) //Considero una alla volta le facce del poliedro originale "precedenti" a f1
			{
				int f2 = polyhedron.Cell2DsId[j];
				vector<int> edges2 = polyhedron.Cell2DsEdges[f2];
				bool found = false; //Variabile booleana indica se ho trovato o no "l'altra faccia" del lato edges1[h]. NOTA: ogni lato è condiviso solo da 2 facce.
				for (int k = 0; k < M; k++)
				{
					if (edges1[h] == edges2[k]){
						dual.Cell1DsId.push_back(newEdge);
						dual.Cell1DsExtrema(0,newEdge) = f1;
						dual.Cell1DsExtrema(1,newEdge) = f2;
						newEdge++; //Incremento id nuovo lato.
						found = true; //Ho trovato faccia in comune.
						//break; //Ho trovato edges1[h] tra i lati di f2 non ha senso considerare gli altri lati di f2.
					}
				}
				if (found) {break;} //Ho trovato la "seconda faccia" di edges1[h], posso passare a condiderare edges1[h+1].
			}
		}
	}
	cout << "newEdge" << newEdge << endl;
	cout << "Cell1DsId: ";
	for (int e : dual.Cell1DsId) {cout << e << " ";};
	cout << endl;
	cout << dual.Cell1DsExtrema << endl;
	
	//Iterando sui vertici del  poliedro originale costruisco le facce del duale
	for (int v : polyhedron.Cell0DsId)
	{
		vector<int> adjacentFaces; //Contenitore dove mi salverò gli id delle facce adiacenti a v nel poliedro originale.
		adjacentFaces.reserve(N);
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
					break;
				}
			}
			if (m >= N) {break;}; //Ci sono al più N facce adiacenti a v.
		}
		dual.Cell2DsId.push_back(v); //v corrisponde ad una faccia del duale.
		dual.Cell2DsVertices.push_back(adjacentFaces); //adjacentFaces corrispondono ai vertici della faccia v nel duale.
		
		
		vector<int> faceEdges; //Contenitore dove vado a salvare gli id dei lati della faccia v (nel duale).
		faceEdges.reserve(N);
		int noE = 0; //Numero di lati di v trovati
		for (auto it1 = adjacentFaces.begin(); it1 != adjacentFaces.end(); it1++) //Scorro i vertici della faccia v
		{
			for (auto it2 = adjacentFaces.begin(); it2 < it1 and it2 != adjacentFaces.end(); it2++)
			{
				Vector2i ex1;
				ex1 << *it1, *it2;
				Vector2i ex2;;
				ex2 << *it2, *it1;
				for (int j = 0; j < E; j++) //Scorro i lati del duale
				{
					if (ex1 == dual.Cell1DsExtrema.col(j) or ex2 == dual.Cell1DsExtrema.col(j))
					{
						faceEdges.push_back(j);
						noE++;
						break;
					};
				}
			}
			if (noE >= N){break;}
		}
		dual.Cell2DsEdges.push_back(faceEdges);
	}
	//
	cout << "Cell2DsId: ";
	for (int e : dual.Cell2DsId) {cout << e << " ";};
	cout << endl;
	cout << "vertici cell2d:" << endl;
	for (vector<int> v : dual.Cell2DsVertices) {
		for (int w : v) {cout << w << " ";};
			cout << endl;
			};
	cout << endl;
	cout << "lati cell2d:" << endl;
	for (vector<int> v : dual.Cell2DsEdges) {
		for (int w : v) {cout << w << " ";};
			cout << endl;
			};
	cout << endl;
	//
	//Aggiorno valori Cell3Ds
	dual.Cell3DsEdges.push_back(dual.Cell1DsId);
	dual.Cell3DsVertices.push_back(dual.Cell1DsId);
	dual.Cell3DsFaces.push_back(dual.Cell2DsId);
return true;	
};


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