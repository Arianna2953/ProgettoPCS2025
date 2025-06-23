#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"
#include "Import.hpp"
#include "Export.hpp"
#include "Triangulation.hpp"
#include "Dual.hpp"
#include "ShortPath.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

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
	cout << "Valori inseriti: p = " << p << ", q = " << q << ", b = " << b << ", c = " << c ;
	
	if(pathToFind){
		cout << ", v0 = " << v0 << ", v1 = " << v1;
	}	
	cout << endl;	
    
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
		cout << "Poliedro regolare di base: Tetraedro (p=3, q=3)" << endl;	
    } 
    else if (p == 3 && q == 4) {
        file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";
		cout << "Poliedro regolare di base: Ottaedro (p=3, q=4)" << endl;
    } 
    else if (p == 3 && q == 5) {
        file0Ds = "../PlatonicSolid/icosahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/icosahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/icosahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/icosahedron/Cell3Ds.txt";
		cout << "Poliedro regolare di base: Icosaedro (p=3, q=5)" << endl;
    } 
    else if (p == 4 && q == 3) {
		//importo il poliedro di base, triangolo, proietto e poi calcolo il duale del poliedro risultante
		
		file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";		
        cout << "Poliedro regolare di base: Ottaedro (p=3, q=4)" << endl;
		dualize = true;
    } 
    else if (p == 5 && q == 3) {
		//importo il poliedro di base, triangolo, proietto e poi calcolo il duale del poliedro risultante
		
		file0Ds = "../PlatonicSolid/icosahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/icosahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/icosahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/icosahedron/Cell3Ds.txt";		
        cout << "Poliedro regolare di base: Icosaedro (p=3, q=5)" << endl;	
		dualize = true;
    } 
    else {
        cerr << "ATTENZIONE - Combinazione di valori di p e q non valida. Nessun poliedro regolare corrispondente." << endl;
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
		if(dualize == true){
			TriangulationTypeI(regularPolyhedron, toDualize,q,p,n);
			cout << "Effettuo la triangolazione di tipo 1" << endl;
			DualConstructor(toDualize,mesh);
			cout << "Effettuo la dualizzazione del poliedro triangolato" << endl;
		}
		else{
			TriangulationTypeI(regularPolyhedron, mesh,p,q,n);
			cout << "Effettuo la triangolazione di tipo 1" << endl;
		}
	}
	else if(b==c && b!=0){
		n = b;
		if(dualize == true){
			TriangulationTypeII(regularPolyhedron, toDualize,n);
			cout << "Effettuo la triangolazione di tipo 2" << endl;
			DualConstructor(toDualize,mesh);
			cout << "Effettuo la dualizzazione del poliedro triangolato" << endl;
		}
		else{
			TriangulationTypeII(regularPolyhedron, mesh,n);
			cout << "Effettuo la triangolazione di tipo 2" << endl;
		}
	}
	else {
		cerr << "ATTENZIONE - Combinazione di valori b e c non valida" << endl;
	}
	
	
	if(pathToFind){
		
		//trovo il percorso minimo sul poliedro tra i vertici v0 e v1
		cout << "CALCOLO DEL CAMMINO MINIMO" << endl;
		cout << "Id vertice di partenza: " << v0 << endl;
		cout << "Id vertice d'arrivo: " << v1 << endl;
		
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