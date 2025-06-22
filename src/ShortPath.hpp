#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
//Data una mesh e gli id di due vertici, restituisce l'id del lato che li collega. Se tale lato non esiste restituisce -1.	
//mesh: PolyhedralMesh struct contenente i dati relativi a vertici, lati e facce di una mesh(di un poliedro).
//v0: id del primo nodo.
//v1: id del secondo nodo.
int FindEdge(const PolyhedralMesh& mesh, const int& v0, const int& v1);

//Data una mesh, costruisce la lista di adiacenza dei suoi nodi.
//mesh: PolyhedralMesh struct contenente i dati relativi a vertici, lati e facce di una mesh(di un poliedro).
//adjList: vettore di liste di numeri interi che rappresenta la lista di adiacenza dei nodi della mesh, adjList[i] è una lista  contenente gli id dei nodi adiacenti al nodo i.
void CreateAdjacencyList(const PolyhedralMesh& mesh, vector<list<int>>& adjList);

//Data una mesh, costruisce la matrice (simmetrica) dei pesi associata assegnando agli elementi (i,j) e (j,i) il valore corrispondente alla distanza tra il nodo i e il nodo j.
//mesh: PolyhedralMesh struct contenente i dati relativi a vertici, lati e facce di una mesh(di un poliedro).
//weights: matrice simmetrica dove vado a salvare in posizione (i,j) il peso associato all'arco i-j, se l'arco i-j non esiste il peso corrispondente avrà valore infinito.
void CreateWeightsMatrix(const PolyhedralMesh& mesh, MatrixXd& weights);

//Dato l'id di un nodo di partenza e l'id di un nodo d'arrivo, visita un grafo pesato secondo l'algoritmo di Dijkstra, 
//calcolando man mano la distanza minima dei nodi visitati dal nodo d'inizio e il loro predecessore nel percorso minimo,
//quando visita il nodo destinazione si ferma.
//adjList: lista di adiacenza dei nodi del grafo.
//s: id sel nodo d'inizio.
//d: id del nodo d'arrivo.
//weights: matrice dei pesi (simmetrica), in posizione (i,j) si trova il peso dell'arco i-j, se l'arco i-j non esiste il peso associato è infinito.
//pred: vettore dei predecessori, perd[i] sarà uguale all'id del nodo che viene prima del nodo i nel percoso minimo, se il predecessore di nodo i è ignodo pred[i] sarà pari a 1.
//dist: vettore delle distanza, dist[i] sarà uguale alla distanza minima del nodo i (dal nodo di partenza), se tale distanza è ignota, dist[i] sarà pari a infinito.
void ComputeDistances(const vector<list<int>>& adjList, const int& s, const int& d, const MatrixXd& weights, vector<int>& pred, vector<double>& dist);

//Data una mesh, l'id di un nodo di partenza e l'id di un nodo d'arrivo, interpreta la mesh come grafo pesato (con la lunghezza dei lati come pesi),
//calcola la lunghezza del camminino minimo tra il nodo di partenza e il nodo di arrivo e il numero di lati che lo compongono.
//Assegna ai nodi e ai lati percorsi dal cammino minimo la proprietà  ShortPath=1.
//mesh: PolyhedralMesh struct contenente i dati relativi a vertici, lati e facce di una mesh(di un poliedro).
//sourceNode: id nodo di partenza.
//destinationNode: id nodo d'arrivo.
//numEdges: numero di lati che compongono il percorso minimo.
//pathLenght: lungheza del percorso minimo.
bool FindShortestPath(PolyhedralMesh& mesh, const int& sourceNode, const int& destinationNode, unsigned int& numEdges, double& pathLenght);

}
