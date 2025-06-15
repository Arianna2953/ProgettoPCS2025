/*
#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
//Triangola le facce del poliedro di partenza, a partire da una divisione dei lati in parti uguali, creando una griglia (?) ????
//polyOld: PolyhedralMesh struct che rappresenta il poliedro di partenza, le cui facce vengono triangolate dalla funzione
//polyNew: PolyhedralMesh struct in cui viene salvato il poliedro con le facce triangolate, e i cui vertici sono proiettati su una sfera di raggio unitario
//p & q: numeri input relativi alla tipologia di poliedro
//n: numero di segmenti in cui viene diviso ogni lato
void TriangulationTypeI(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& p, const int& q, const int& n); 

//Triangola le facce del poliedro di partenza, a partire da una divisione dei lati in parti uguali, creando una griglia che considera anche l'unione con i baricentri dei nuovi triangoli creati(?) ????
//(da capire bene la logica prima di commentare...)
//polyOld: PolyhedralMesh struct che rappresenta il poliedro di partenza, le cui facce vengono triangolate dalla funzione
//polyNew: PolyhedralMesh struct in cui viene salvato il poliedro con le facce triangolate, e i cui vertici sono proiettati su una sfera di raggio unitario
//n: numero di segmenti in cui viene diviso ogni lato
void TriangulationTypeII(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& n); 

/*Durante la creazione della griglia triangolare, controlla se un lato è già stato inserito nella struttura: il lato considerato viene inserito se non ancora presente
  e la funzione restituisce il valore del suo id, se invece è già stato inserito, viene restituito il valore dell'id del lato corrispondente già presente nella mesh*/
//poly: PolyhedralMesh struct in cui si vuole valutare l'esistenza e l'eventuale inserimento del lato edge
//edge: lato di cui si vuole valutare l'inserimento nella PolyhedralMesh
//id_edge: id del lato che viene aggiornato ogni volta che viene aggiunto un nuovo lato alla struttura
int CheckAddEdges(PolyhedralMesh& poly, const Vector2i& edge, int& id_edge);

//Controlla se un nuovo punto che si vuole inserire sia già stato inserito tra i vertici del poliedro considerato
//poly: PolyhedralMesh struct in cui si vuole valutare l'esistenza e l'eventuale inserimento del lato edge
//vertex: vertice di cui si vuole valutare l'inserimento nella PolyhedralMesh
//id_vert: id del vertice che viene aggiornato ogni volta che viene aggiunto un nuovo lato alla struttura
int CheckAddVertices(PolyhedralMesh& poly, const Vector3d& vertex, int& id_vert);

}
*/