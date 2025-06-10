#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"
#include "ImportExport.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

/*
---> per shortpath Ele
//creo una lista di adiacenza nella forma di un vettore di liste
vector<list<int>> CreateAdjacencyList(const PolyhedralMesh& mesh) {
    vector<list<int>> adjacencyList(mesh.NumCell0Ds);

    for (int i = 0; i < mesh.NumCell1Ds; i++){
        int idFrom = mesh.Cell1DsExtrema(0, i);
        int idTo = mesh.Cell1DsExtrema(1, i);

        adjacencyList[idFrom].push_back(idTo);
        adjacencyList[idTo].push_back(idFrom);
    }
	for (size_t i = 0; i < adjacencyList.size(); ++i) {
        cout << "Vertice " << i << " → ";
        for (int neighbor : adjacencyList[i]) {
            cout << neighbor << " ";
        }
	cout << endl;}
    return adjacencyList;
}*/

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
    int v0 = 0, v1 = 0;
    if (argc == 7) {
        v0 = stoi(argv[5]);
        v1 = stoi(argv[6]);
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
		//importo il poliedro di base, triangolo, proietto e poi calcolo il duale??
		
		file0Ds = "../PlatonicSolid/octahedron/Cell0Ds.txt";
		file1Ds = "../PlatonicSolid/octahedron/Cell1Ds.txt";
		file2Ds = "../PlatonicSolid/octahedron/Cell2Ds.txt";
		file3Ds = "../PlatonicSolid/octahedron/Cell3Ds.txt";		
        cout << "Poliedro regolare di base: Ottaedro - bisogna poi farne il duale" << endl;
		dualize = true;
    } 
    else if (p == 5 && q == 3) {
		//importo il poliedro di base, triangolo, proietto e poi calcolo il duale??
		
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
			cerr << "file non trovato" << endl;
			return 1;
		} 
	
	int n;
	
	//imposto i vari casi in base al valore di b e c
	if((b==0 && c >=1) || (b>=1 && c==0)){
		n = max(b,c);
		cout << "Triangolazione di 'tipo 1'" << endl;
		if(dualize == true){
			TriangulationTypeI(regularPolyhedron, toDualize,p,q,n);
			DualConstructor(toDualize,mesh);
		}
		else{
			TriangulationTypeI(regularPolyhedron, mesh,p,q,n);
		}
	}
	else if(b==c && b!=0){
		n = b;
		//triangolazione tipo 2
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
	

	/*{
    vector<Gedim::UCDProperty<double>> cell2Ds_properties(1);

    cell2Ds_properties[0].UnitLabel = "-";
    cell2Ds_properties[0].NumComponents = 1;

    // Crea lista dei vertici per ogni poligono
    std::vector<std::vector<unsigned int>> polygons_vertices(mesh.NumCell2Ds);
    for (int i = 0; i < mesh.NumCell2Ds; ++i)
    {
        for (int v : mesh.Cell2DsVertices[i])
            polygons_vertices[i].push_back(static_cast<unsigned int>(v));
    }

    // Crea matrice punti 3D: (3 x NumCell0Ds)
    Eigen::MatrixXd points3D(3, mesh.NumCell0Ds);
    for (int i = 0; i < mesh.NumCell0Ds; ++i)
    {
        points3D(0, i) = mesh.Cell0DsCoordinates(0, i);
        points3D(1, i) = mesh.Cell0DsCoordinates(1, i);
        points3D(2, i) = mesh.Cell0DsCoordinates(2, i);
    }

    // Materiali (opzionale, azzerati qui)
    Eigen::VectorXi materials = Eigen::VectorXi::Zero(mesh.NumCell2Ds);

    utilities.ExportPolygons("./Cell2Ds.inp",
                             points3D,
                             polygons_vertices,
                             {},
                             cell2Ds_properties,
                             materials);
}*/
	/*//visualizzazione a terminale
    cout << "=== Cell0Ds (Vertices) ===\n";
    cout << "NumCell0Ds: " << mesh.NumCell0Ds << "\n";
    cout << "Cell0DsId: ";
    for (int id : mesh.Cell0DsId) cout << id << " ";
    cout << "\nCoordinates:\n" << mesh.Cell0DsCoordinates << "\n";
    cout << "ShortPathCell0Ds:\n";
    for (const auto& [shortpath, ids] : mesh.ShortPathCell0Ds) {
        cout << "  ShortPath " << shortpath << ": ";
        for (int id : ids) cout << id << " ";
        cout << "\n";
    }

    cout << "\n=== Cell1Ds (Edges) ===\n";
    cout << "NumCell1Ds: " << mesh.NumCell1Ds << "\n";
    cout << "Cell1DsId: ";
    for (int id : mesh.Cell1DsId) cout << id << " ";
    cout << "\nExtrema:\n" << mesh.Cell1DsExtrema << "\n";
    cout << "ShortPathCell1Ds:\n";
    for (const auto& [shortpath, ids] : mesh.ShortPathCell1Ds) {
        cout << "  ShortPath " << shortpath << ": ";
        for (int id : ids) cout << id << " ";
        cout << "\n";
    }

    cout << "\n=== Cell2Ds (Faces) ===\n";
    cout << "NumCell2Ds: " << mesh.NumCell2Ds << "\n";
    cout << "Cell2DsId: ";
    for (int id : mesh.Cell2DsId) cout << id << " ";
    cout << "\nVertices per Cell2D:\n";
    for (size_t i = 0; i < mesh.Cell2DsVertices.size(); ++i) {
        cout << "  Face " << i << ": ";
        for (int v : mesh.Cell2DsVertices[i]) cout << v << " ";
        cout << "\n";
    }
    cout << "Edges per Cell2D:\n";
    for (size_t i = 0; i < mesh.Cell2DsEdges.size(); ++i) {
        cout << "  Face " << i << ": ";
        for (int e : mesh.Cell2DsEdges[i]) cout << e << " ";
        cout << "\n";
    }

    cout << "\n=== Cell3Ds (Volumes) ===\n";
    cout << "NumCell3Ds: " << mesh.NumCell3Ds << "\n";
    cout << "Cell3DsId: ";
    for (int id : mesh.Cell3DsId) cout << id << " ";
    cout << "\nVertices per Cell3D:\n";
    for (size_t i = 0; i < mesh.Cell3DsVertices.size(); ++i) {
        cout << "  Volume " << i << ": ";
        for (int v : mesh.Cell3DsVertices[i]) cout << v << " ";
        cout << "\n";
    }
    cout << "Edges per Cell3D:\n";
    for (size_t i = 0; i < mesh.Cell3DsEdges.size(); ++i) {
        cout << "  Volume " << i << ": ";
        for (int e : mesh.Cell3DsEdges[i]) cout << e << " ";
        cout << "\n";
    }
    cout << "Faces per Cell3D:\n";
    for (size_t i = 0; i < mesh.Cell3DsFaces.size(); ++i) {
        cout << "  Volume " << i << ": ";
        for (int f : mesh.Cell3DsFaces[i]) cout << f << " ";
        cout << "\n";
    }*/


}