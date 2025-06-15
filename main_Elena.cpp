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

int FindEdge(const PolyhedralMesh& mesh, const int& v0, const int& v1){
	int edgeId = -1;
	for (int j = 0; j < mesh.NumCell1Ds; j++){
		int u0 = mesh.Cell1DsExtrema(0,j);
		int u1 = mesh.Cell1DsExtrema(1,j);
		if ( (v0 == u0 && v1 == u1) || (v0 == u1 && v1 == u0) ){
			edgeId = j;
			break;
		}
	}
	return edgeId;
}

void CreateAdjacencyList(const PolyhedralMesh& mesh, vector<list<int>>& adjList){
    for (int i = 0; i < mesh.NumCell1Ds; i++){
        int idFrom = mesh.Cell1DsExtrema(0, i);
        int idTo = mesh.Cell1DsExtrema(1, i);

        adjList[idFrom].push_back(idTo);
        adjList[idTo].push_back(idFrom);
    }
    return;
}

void CreateWheightsMatrix(const PolyhedralMesh& mesh, MatrixXd& weights){
	for (int e = 0; e < mesh.NumCell1Ds; e++){
		int idV1 = mesh.Cell1DsExtrema(0,e);
		int idV2 = mesh.Cell1DsExtrema(1,e);
		Vector3d coordV1 = mesh.Cell0DsCoordinates.col(idV1);
		Vector3d coordV2 = mesh.Cell0DsCoordinates.col(idV2);
		double edgeLenght = (coordV1-coordV2).norm();
		weights(idV1, idV2) = edgeLenght;
		weights(idV2, idV1) = edgeLenght;
		}
	return;
}

void ComputeDistances(const vector<list<int>>& adjList, const int& s, const int& d, const MatrixXd& weights, const int& V, vector<int>& pred, vector<double>& dist){
	//Riempio il vettore dei predecessori e quello delle distanze.
	for (int i = 0; i < V; i++){
		pred[i] = -1;
		dist[i] = INFINITY;	
	}
	
	//Inizializzo il predecessore e la distanza del nodo sorgente s.
	pred[s] = s;
	dist[s] = 0;
	
	//Creo un vettore di n elementi per tenere traccia dei nodi  già  visitati.
	vector<int> visitedNodes(V,0); //inizializzo a zero
	//for (int i = 0; i < V; i++) {visitedNodes.push_back(0);} 
	
	//Implemento la priority queue con un heap partendo da un vettore di tuple (distanza_nodo_i, nodo_i).
	vector<pair<double,int>> priorityQueue;
	for (int i = 0; i < V; i++){
		priorityQueue.push_back(pair(dist[i], i));
	}
	make_heap(priorityQueue.begin(), priorityQueue.end(), greater<>{});
		
	for(int i = 0; i < V-1; i++){ //NOTA: uso un for per i che va da 0 a V-2 perché posso visitare al massimo V nodi e visitare l'ultimo (il più lontano) sarebbe superfluo
		const int u = get<1>(priorityQueue.front()); //Leggo l'id del nodo con distanza minima.
		
		if (u == d) {break;} //Controllo se ho raggiunto il nodo destinazione, se sì esco  dal ciclo for.
		
		visitedNodes[u]  = 1; //Aggiorno visitedNodes, segno che il nodo u è stato visitato.
		
		//Dequeue
		pop_heap(priorityQueue.begin(), priorityQueue.end());
		priorityQueue.pop_back();
		
		for (int w : adjList[u]){
			if (visitedNodes[w] != 0) {continue;} //Se ho già visitato il nodo w, passo al sucessivo
			if (dist[w] > dist[u] + weights(u,w)){
				pred[w] = u;
				dist[w] = dist[u] + weights(u,w);
				priorityQueue.push_back(pair(dist[w],w));
				push_heap(priorityQueue.begin(), priorityQueue.end(), greater<>{});
			}
		}
		
		cout << "Nodo estratto:" << u << endl;
		for (int i = 0; i < V; i++) {cout << visitedNodes[i] << "  ";}
		cout << endl;
		for (int i = 0; i < V; i++) {cout << dist[i] << "  ";}
		cout << endl;
		for (int i = 0; i < V; i++) {cout << pred[i] << "  ";}
		cout << endl;
	
	}
	 
	return;
}

bool FindShortestPath(PolyhedralMesh& mesh, const int& sourceNode, const int& destinationNode){
	for ( int i : mesh.Cell0DsId) {cout << i << " ";}
	cout << endl;
	
	const int V = mesh.NumCell0Ds;
	
	if (sourceNode < 0 || sourceNode >= V || destinationNode < 0 || destinationNode >= V ){
		cerr << "Id vertici non validi." << endl;
		return false;
	}
	
	cout << "Id vertice di partenza: " << sourceNode << endl;
	cout << "Id vertice d'arrivo: " << destinationNode << endl;
	
	//Creo la lista di adiacenza dei vertici del poliedro come vettore di liste.
	vector<list<int>> adjacencyList(V);
	CreateAdjacencyList(mesh, adjacencyList);
	
	//Creo la matrice dei pesi.
	MatrixXd weightsEdges = MatrixXd::Ones(V,V)*INFINITY;
	CreateWheightsMatrix(mesh, weightsEdges);
	
	vector<int> predecessors(V);
	vector<double> distances(V);
	ComputeDistances(adjacencyList, sourceNode, destinationNode, weightsEdges, V, predecessors, distances);
	
	for (int v : predecessors) {cout << v << " ";}
	cout << endl;
	
	//Calcolo il numero di lati nel percorso minimo 
	//e aggiorno la proprietà ShortPath dei nodi e dei lati che compongono il percorso.
	int ShortPath = 1;
	
	unsigned int numEdges =0;
	int currentNode = destinationNode;
	
	// Ricostruisco il cammino minimo andando indietro
	while (currentNode != sourceNode) {
		mesh.ShortPathCell0Ds[ShortPath].push_back(currentNode); 
		
		int prevNode = predecessors[currentNode]; 
		int edgeId = FindEdge(mesh, currentNode, prevNode);
		//
		cout << "v: " << currentNode << " prev: " << prevNode << " lato: " << edgeId << endl;
		//
		if (edgeId < 0 || edgeId >= mesh.NumCell1Ds) {
			cerr << "Id lato non valido trovato." << endl;
			return false;
		}
		mesh.ShortPathCell1Ds[ShortPath].push_back(edgeId);
		
		currentNode = prevNode;
		numEdges++;
	}
	// Aggiungo anche il nodo sorgente alla lista dei nodi nel cammino minimo
	mesh.ShortPathCell0Ds[ShortPath].push_back(sourceNode);
	
	//rimuovo i vertici e i lati dal vettore associato a ShortPath = 0
	for (int v : mesh.ShortPathCell0Ds[ShortPath]) {
			mesh.ShortPathCell0Ds[0].remove(v);  // rimuove dalla lista "default" con shortpath = 0
		}
	
	for (int e : mesh.ShortPathCell1Ds[ShortPath]) {
			mesh.ShortPathCell1Ds[0].remove(e);  // rimuove dalla lista "default" con shortpath = 0
		}
	
	// Inverto l'ordine dei nodi e spigoli (perché sono salvati "al contrario")
	//Voglio il cammino nel "verso giusto": da sourceNode a destinationNode
	reverse(mesh.ShortPathCell0Ds[ShortPath].begin(), mesh.ShortPathCell0Ds[ShortPath].end());
	reverse(mesh.ShortPathCell1Ds[ShortPath].begin(), mesh.ShortPathCell1Ds[ShortPath].end());
	
	//per vedere quali nodi e archi sono parte del percorso minimo (check a terminale... si può poi togliere)
	cout << "Cammino minimo (nodi): ";
	for (int v : mesh.ShortPathCell0Ds[1]) {cout << v << " ";}
	cout << endl;

	cout << "Cammino minimo (spigoli): ";
	for (int e : mesh.ShortPathCell1Ds[1]) {cout << e << " ";}
	cout << endl;
	
	cout << "Numero di lati percorso minimo: " << numEdges << endl;
	cout << "Lunghezza percorso minimo: " << distances[destinationNode] << endl;
	
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
	
	//trovo il percorso minimo sul poliedro tra i vertici v0 e v1 
	if(v0 >= 0 && v1 >= 0){
		FindShortestPath(mesh,v0,v1);	
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