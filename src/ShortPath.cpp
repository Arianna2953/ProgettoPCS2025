#include "ShortPath.hpp"
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
#include <queue>

namespace PolyhedralLibrary
{

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

void CreateAdjacencyList(const PolyhedralMesh& mesh, vector<list<int>>& adjList){
    for (int i = 0; i < mesh.NumCell1Ds; i++){
        int idFrom = mesh.Cell1DsExtrema(0, i);
        int idTo = mesh.Cell1DsExtrema(1, i);

        adjList[idFrom].push_back(idTo);
        adjList[idTo].push_back(idFrom);
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

void CreateWeightsMatrix(const PolyhedralMesh& mesh, MatrixXd& weights){
	for (int j = 0; j < mesh.NumCell1Ds; j++){
		int idV1 = mesh.Cell1DsExtrema(0,j);
		int idV2 = mesh.Cell1DsExtrema(1,j);
		Vector3d coordV1 = mesh.Cell0DsCoordinates.col(idV1);
		Vector3d coordV2 = mesh.Cell0DsCoordinates.col(idV2);
		double edgeLenght = (coordV1-coordV2).norm();
		weights(idV1, idV2) = edgeLenght;
		weights(idV2, idV1) = edgeLenght;
		}
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

void ComputeDistances(const vector<list<int>>& adjList, const int& s, const int& d, const MatrixXd& weights, vector<int>& pred, vector<double>& dist){
	
	const int V = adjList.size(); //Numero di nodi del grafo.
	
	//Riempio il vettore dei predecessori e quello delle distanze.
	pred.assign(V,-1);
	dist.assign(V,INFINITY);
	
	//Inizializzo il predecessore e la distanza del nodo sorgente s.
	pred[s] = s;
	dist[s] = 0;
	
	//Creo un vettore di n elementi per tenere traccia dei nodi visitati (0: non visitato; 1: visitato).
	vector<int> visited(V,0);
	
	//Priority queue. Min heap : (distanza, nodo)
	priority_queue<pair<double,int>,vector<pair<double,int>>,greater<>> pq;
	for (int i = 0; i < V; i++){
		pq.push(pair(dist[i],i));
	}
	
	while(!pq.empty()){ 
		const int u = get<1>(pq.top()); //Leggo l'id del nodo con distanza minima.
		
		pq.pop();//Dequeue, rimuovo il nodo u dalla priority queue.
		
		if (u == d) {break;} //Controllo se ho raggiunto il nodo destinazione, se sì esco dal ciclo for.
		if (visited[u] != 0) {continue;} //Se ho già visitato il nodo u passo all'iterazione successiva
		
		visited[u] = 1; //Aggiorno visited, segno che il nodo u è stato visitato.
		
		for (int w : adjList[u]){
			if (visited[w] != 0) {continue;} //Se ho già visitato il nodo w, passo al sucessivo.
			if (dist[w] > dist[u] + weights(u,w)){
				pred[w] = u;
				dist[w] = dist[u] + weights(u,w);
				pq.push(pair(dist[w],w));
			}
		}
	}
	 
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

bool FindShortestPath(PolyhedralMesh& mesh, const int& sourceNode, const int& destinationNode, unsigned int& numEdges, double& pathLenght){
	
	/*cout << "CALCOLO DEL CAMMINO MINIMO" << endl;
	cout << "Id vertice di partenza: " << sourceNode << endl;
	cout << "Id vertice d'arrivo: " << destinationNode << endl;*/
	
	const int V = mesh.NumCell0Ds; //Numero di nodi del grafo.
	
	//Controllo che gli id dei nodi di partenza e destinazione siano validi.
	if (sourceNode < 0 || sourceNode >= V || destinationNode < 0 || destinationNode >= V ){
		cerr << "ATTENZIONE: Id vertici inseriti non validi." << endl;
		return false;
	}
	
	//Creo la lista di adiacenza dei vertici del poliedro come vettore di liste di interi.
	vector<list<int>> adjacencyList(V);
	CreateAdjacencyList(mesh, adjacencyList);
	
	//Creo la matrice dei pesi.
	MatrixXd weightsEdges = MatrixXd::Ones(V,V)*INFINITY;
	CreateWeightsMatrix(mesh, weightsEdges);
	
	//Calcolo la distanza minima con l'algoritmo di Dijkstra.
	vector<int> predecessors;
	vector<double> distances;
	ComputeDistances(adjacencyList, sourceNode, destinationNode, weightsEdges, predecessors, distances);
	
	pathLenght = distances[destinationNode]; 
	
	//Calcolo il numero di lati nel percorso minimo e aggiorno la proprietà ShortPath dei nodi e dei lati che compongono il percorso.
	int ShortPath = 1;
	int currentNode = destinationNode;
	
	// Parto dal nodo destinazione e ricostruisco il cammino minimo andando indietro.
	while (currentNode != sourceNode) {
		mesh.ShortPathCell0Ds[ShortPath].push_back(currentNode); 
		
		int prevNode = predecessors[currentNode]; 
		int edgeId = FindEdge(mesh, currentNode, prevNode);
		if (edgeId < 0 || edgeId >= mesh.NumCell1Ds) {
			cerr << "Id lato non valido trovato." << endl;
			return false;
		}
		mesh.ShortPathCell1Ds[ShortPath].push_back(edgeId);
		
		currentNode = prevNode;
		numEdges++;
	}
	// Aggiorno la proprietà  ShortPath del nodo di partenza.
	mesh.ShortPathCell0Ds[ShortPath].push_back(sourceNode);
	
	//Rimuovo i vertici e i lati del cammino minimo dalla lista "di default" associata a ShortPath = 0.
	for (int v : mesh.ShortPathCell0Ds[ShortPath]) {
			mesh.ShortPathCell0Ds[0].remove(v);  
		}
	for (int e : mesh.ShortPathCell1Ds[ShortPath]) {
			mesh.ShortPathCell1Ds[0].remove(e);
		}
	
	
	//Inverto l'ordine dei nodi e spigoli (perché sono salvati "al contrario")
	//Voglio il cammino nel "verso giusto": da sourceNode a destinationNode
	reverse(mesh.ShortPathCell0Ds[ShortPath].begin(), mesh.ShortPathCell0Ds[ShortPath].end());
	reverse(mesh.ShortPathCell1Ds[ShortPath].begin(), mesh.ShortPathCell1Ds[ShortPath].end());
	
	/*//Stampo a terminale i nodi e i lati del percorso minimo
	cout << "Cammino minimo (nodi): ";
	for (int v : mesh.ShortPathCell0Ds[1]) {cout << v << " ";}
	cout << endl;
	cout << "Cammino minimo (spigoli): ";
	for (int e : mesh.ShortPathCell1Ds[1]) {cout << e << " ";}
	cout << endl;
	*/
	
	return true;
}


}

