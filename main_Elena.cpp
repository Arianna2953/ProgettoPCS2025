#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"

//#include "ImportExport.hpp"
//#include "Utils.hpp"
#include "Import.hpp"
#include "Export.hpp"
#include "Triangulation.hpp"
//#include "Dual.hpp"
#include "ShortPath.hpp"

#include "UCDUtilities.hpp"

#include <queue>
#include <cmath>
#include <list>
#include <map>
#include <vector>
#include <set>

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;


bool DualConstructor(const PolyhedralMesh& polyhedron, PolyhedralMesh& dual){	
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
	for (int id : polyhedron.Cell2DsId){
		Vector3d barCoords; //Contenitore dove salverò le coordinate del baricentro
		
		//Ricavo id e coordinate dei vertici della faccia
		const vector<int>& faceVertices = polyhedron.Cell2DsVertices[id];
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
		dual.Cell0DsId.push_back(id);
		dual.Cell0DsCoordinates.col(id) << barCoords;
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
			//if (m >= Ndual) {break;}; //Ci sono al più Ndual facce adiacenti a v.
		}
		vector<int>& faceVerticesDual = adjacentFaces; //Quelli che nel poliedro originale erano gli id delle facce adiacenti al vertice v, nel duale sono gli id dei vertici della faccia con id == v.
		
		
		vector<int> faceEdgesDual(m); 
		//faceEdgesDual.reserve(Ndual);
		vector<int> fVDSorted(m);
		//fVDSorted.reserve(Ndual);
		
		bool isClosed = false;
		int u0 = faceVerticesDual[0];
		int u1 = faceVerticesDual[0];
		
		for(int i = 0; i < m; i++){
			bool foundE = false;
			for (int u2 : faceVerticesDual){
				if (u2 == u1 || u2 == u0) {continue;}
				for(int j = 0; j < E; j++){
					int ex1 = dual.Cell1DsExtrema(0,j);
					int ex2 = dual.Cell1DsExtrema(1,j);
					if((u1 == ex1 && u2 == ex2)||(u1 == ex2 && u2 == ex1)){
						faceEdgesDual.push_back(j);
						fVDSorted.push_back(u1);
						if(u2 == fVDSorted[0]){isClosed = true;}
						u0 = u1;
						u1 = u2;
						foundE = true;
						break;
					}
				}
				if (foundE) {break;}
			}
		}
		if(!isClosed) {
			cerr << "ATTENZIONE: faccia del duale aperta." << endl;
			return false;
		}
		
		dual.Cell2DsId.push_back(v); //ATTENZIONE: associo alla nuova faccia del duale lo stesso id del vertice del poliedro di partenza corrispondente
		dual.Cell2DsVertices.push_back(fVDSorted); 
		dual.Cell2DsEdges.push_back(faceEdgesDual);
		
	}
	
	
	//Aggiorno valori Cell3Ds
	dual.Cell3DsVertices.push_back(dual.Cell0DsId);
	dual.Cell3DsEdges.push_back(dual.Cell1DsId);
	dual.Cell3DsFaces.push_back(dual.Cell2DsId);
	return true;	
}


int main(int argc, char* argv[]) {
    // Verifica che ci siano 5 o 7 argomenti (5 se non ci sono v0 e v1, 7 se ci sono v0 e v1)
    if (argc != 5 && argc != 7) {
        cerr << "Indicare i parametri nella forma: ./programma p q b c [v0 v1]\n";
        return 1;
    }

    // Cicla su tutti gli argomenti (salta argv[0] che è il nome del programma)
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        // Controlla se l'argomento è un numero intero positivo
        if (arg.empty() || arg[0] == '-' || !std::all_of(arg.begin(), arg.end(), ::isdigit)) {
            cerr << "Errore: argomento " << i << " (" << arg << ") non è un numero intero positivo valido\n";
            return 1;
        }
    }

    // Se gli argomenti sono validi, convertili in numeri interi
    int p = stoi(argv[1]);
    int q = stoi(argv[2]);
    int b = stoi(argv[3]);
    int c = stoi(argv[4]);
    int v0 = -1, v1 = -1;
	bool pathToFind = false;
    if (argc == 7) {
        v0 = stoi(argv[5]);
        v1 = stoi(argv[6]);
		pathToFind = true;
    }

    cout << "Valori inseriti: p = " << p << ", q = " << q << ", b = " << b << ", c = " << c << ", v0 = " << v0 << ", v1 = " << v1 << endl;
    
    string file0Ds;
    string file1Ds;
    string file2Ds;
    string file3Ds;
	
	PolyhedralMesh regularPolyhedron;
	PolyhedralMesh toDualize;
	bool dualize = false;
	PolyhedralMesh mesh;

    // Controllo combinazioni specifiche di p e q
    if (p == 3 && q == 3) {
		file0Ds = "../PlatonicSolid/tetrahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/tetrahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/tetrahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/tetrahedron/Cell3Ds.txt";
		cout << "Poliedro regolare di base: Tetraedro" << endl;	
    } 
    else if (p == 3 && q == 4) {
        file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
		cout << "Poliedro regolare di base: Ottaedro" << endl;
    } 
    else if (p == 3 && q == 5) {
        file0Ds = "../PlatonicSolid/icosahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/icosahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/icosahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/icosahedron/Cell3Ds.txt";
		cout << "Poliedro regolare di base: Icosaedro" << endl;
    } 
    else if (p == 4 && q == 3) {
		//importo il poliedro di base, triangolo, proietto e poi calcolo il duale del poliedro risultante
		
		file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";		
        cout << "Poliedro regolare di base: Ottaedro - bisogna poi farne il duale" << endl;
		dualize = true;
    } 
    else if (p == 5 && q == 3) {
		//importo il poliedro di base, triangolo, proietto e poi calcolo il duale del poliedro risultante
		
		file0Ds = "../PlatonicSolid/icosahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/icosahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/icosahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/icosahedron/Cell3Ds.txt";		
        cout << "Poliedro regolare di base: Icosaedro - bisogna poi farne il duale" << endl;	
		dualize = true;
    } 
    else {
        cout << "Combinazione di valori di p e q non corrispondente a nessun poliedro regolare." << endl;
		return 1;
    }

	if(!ImportMesh(regularPolyhedron,file0Ds,file1Ds,file2Ds,file3Ds))
		{
			cerr << "File non trovato." << endl;
			return 1;
		} 
	
	int n;
	
	//imposto i vari casi in base al valore di b e c
	if((b==0 && c >=1) || (b>=1 && c==0)){
		n = max(b,c);
		cout << "Triangolazione di 'tipo 1'" << endl;
		if(dualize == true){
			TriangulationTypeI(regularPolyhedron, toDualize,q,p,n);
			DualConstructor(toDualize,mesh);
		}
		else{
			TriangulationTypeI(regularPolyhedron, mesh,p,q,n);
		}
	}
	else if(b==c && b!=0){
		n = b;
		cout << "Triangolazione di 'tipo 2'" << endl;
		if(dualize == true){
			TriangulationTypeII(regularPolyhedron, toDualize,n);
			DualConstructor(toDualize,mesh);
		}
		else{
			TriangulationTypeII(regularPolyhedron, mesh,n);
		}
	}
	else {
		cout << "valori di b e c non validi" << endl;
	}
	
	//trovo il percorso minimo sul poliedro tra i vertici v0 e v1 
	if(pathToFind){
		unsigned int numEdges = 0;
		double pathLenght = 0.0;
		if(!FindShortestPath(mesh,v0,v1,numEdges,pathLenght)){
			cerr << "Impossibile trovare il percorso minimo." << endl;
			return 1;
		}	
		cout << "Numero di lati percorso minimo: " << numEdges << endl;
		cout << "Lunghezza percorso minimo: " << pathLenght << endl;
	}
	
	ExportPolyhedron(mesh);
	
	
//esporto file per visualizzazione
Gedim::UCDUtilities utilities;
    {	vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
        for(const auto &m : mesh.ShortPathCell0Ds)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();

        utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates,
                               cell0Ds_properties);
    }

    {	vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
        for(const auto &m : mesh.ShortPathCell1Ds)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
    }
}

