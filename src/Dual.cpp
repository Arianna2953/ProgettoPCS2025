#include "Dual.hpp"
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

namespace PolyhedralLibrary
{
	
bool DualConstructor(const PolyhedralMesh& polyhedron, PolyhedralMesh& dual){	
	const int V = polyhedron.NumCell2Ds; //Numero facce poliedro originale = numero vertici duale
	const int E = polyhedron.NumCell1Ds; //Numero lati poliedro originale = numero lati duale
	const int F = polyhedron.NumCell0Ds; //Numero vertici poliedro originale = numero facce duale
	
	//Alloco memoria a sufficienza per memorizzare tutti gli elementi del duale.
	dual.NumCell0Ds = V;
	dual.Cell0DsId.reserve(V); 
	dual.Cell0DsCoordinates = MatrixXd::Zero(3,V); 
	
	dual.NumCell1Ds = E;
	dual.Cell1DsId.reserve(E);
	dual.Cell1DsExtrema = MatrixXi::Zero(2,E); 
	
	dual.NumCell2Ds = F;
	dual.Cell2DsId.reserve(F); 
	dual.Cell2DsVertices.reserve(F);
	dual.Cell2DsEdges.reserve(F);
	
	dual.NumCell3Ds = 1;
	dual.Cell3DsId.push_back(0);
	dual.Cell3DsEdges.reserve(1);
	dual.Cell3DsVertices.reserve(1);
	dual.Cell3DsFaces.reserve(1);
	
	const int Npoly = 3; //Numero di vertici/lati per faccia nel poliedro originale 
	const int Ndual = E*2/F; //Numero di vertici/lati per faccia nel poliedro duale (circa)
	
	//Iterando sulle facce del poliedro originale trovo i baricentri, cioè i vertici del poliedro duale
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
	
	//Trovo i lati del duale. 
	//Scorro le facce (f1) del poliedro orginale e ne considero un lato alla volta cercando la faccia (f2) 
	//con cui è condiviso. I baricentri di f1 e f2, nel duale, saranno gli estremi di un lato, che vado ad aggiungere all'elenco.
	int newEdge = 0; //Id lato da aggiungere
	for (int i = 0; i < polyhedron.NumCell2Ds; i++) //Itero sulle facce del poliedro originale
	{
		int f1 = polyhedron.Cell2DsId[i];
		vector<int> edgesf1 = polyhedron.Cell2DsEdges[f1];
		
		for (int h = 0; h < Npoly; h++) //Considero uno alla volta i lati della faccia f1 
		{
			bool found = false; //Variabile booleana indica se ho trovato o no "l'altra faccia" a cui appartiene il lato edgesf1[h]. NOTA: ogni lato è condiviso solo da 2 facce.
			for (int j = 0; j <i; j ++) //Considero una alla volta le facce del poliedro originale "precedenti" a f1
			{
				int f2 = polyhedron.Cell2DsId[j];
				vector<int> edgesf2 = polyhedron.Cell2DsEdges[f2];
				for (int k = 0; k < Npoly; k++) //Itero sui lati di f2
				{
					if (edgesf1[h] == edgesf2[k]){ //Se le facce f1 e f2 hanno un lato in comune nel poliedro originale, aggiungo il lato di estremi f1 e f2 (intesi come baricentri) all'elenco dei lati del duale
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
	
	//Iterando sui vertici del poliedro originale costruisco le facce del duale, ad ogni vertice corrisponderà una faccia che condividerà lo stesso id.
	//Considero un vertice del poliedro originale alla volta e ne trovo le facce adiacenti.
	//Passando al duale, interpreto gli id delle facce trovate come id dei loro baricentri, ovvero come id di vertici del duale.
	//Partendo dal primo vertice trovato, percorro il bordo della faccia, risalendo così agli insiemi ordinati dei suoi vertici e lati. 
	for (int v : polyhedron.Cell0DsId) {
		vector<int> adjacentFaces; //Contenitore in cui salverò gli id delle facce adiacenti al vertice v
		adjacentFaces.reserve(Ndual);
		
		// Trova tutte le facce adiacenti al vertice v
		for (int f = 0; f < polyhedron.NumCell2Ds; ++f) {
			for (int vf : polyhedron.Cell2DsVertices[f]) {
				if (v == vf) {
					adjacentFaces.push_back(f);
					break;
				}
			}
		}

		//Reinterpreto gli id delle facce adiacenti al vertice v nel poliedro originale come id dei vertici della faccia corrispondente nel duale
		const vector<int>& dualFaceUnorderedVertices = adjacentFaces; 

		// Ricostruzione ordinata dei vertici della faccia duale e dei sui lati
		vector<int> orderedVertices;
		vector<int> orderedEdges;
		vector<bool> visited(polyhedron.NumCell2Ds, false); //Segno i vertici già inseriti nel vettore ordinato, per non ripetere

		//Parto dal primo vertice (baricentro della prima faccia) trovato
		int current = dualFaceUnorderedVertices[0];
		visited[current] = true;
		orderedVertices.push_back(current);
		
		//Cerco il vertice successivo adiacente a quello corrente tramite il lato che li unisce (nel duale)
		while (orderedVertices.size() < dualFaceUnorderedVertices.size()) {
			bool found = false;
			for (int u : dualFaceUnorderedVertices) {
				if (visited[u]) continue; //Se ho già visitato u passo al vertice successivo
				
				for (int j = 0; j < E; ++j) { //Ciclo for sui lati del duale
					Vector2i e = dual.Cell1DsExtrema.col(j);
					if ((e[0] == current && e[1] == u) || (e[1] == current && e[0] == u)) { //Il j-esimo lato del duale unisce i vertici current e u? Se sì...
						found = true; //Ho trovato il vertice successivo a current
						visited[u] = true;
						orderedVertices.push_back(u);
						orderedEdges.push_back(dual.Cell1DsId[j]);
						current = u; //Adesso il vertice corrente è u
						break;
					}
				}
				if (found) break;
			}

			if (!found) {
				cerr << "ATTENZIONE: ciclo non chiuso correttamente la faccia " << v <<  " del poliedro duale." << endl;
				break;
			}
		}
		//Aggiungo l'ultimo lato per chiudere la faccia
		if (orderedVertices.size() == dualFaceUnorderedVertices.size()) {
			for (int j = 0; j < E; ++j) {
				Vector2i e = dual.Cell1DsExtrema.col(j);
				if ((e[0] == current && e[1] == orderedVertices[0]) || (e[1] == current && e[0] == orderedVertices[0])) {
					orderedEdges.push_back(dual.Cell1DsId[j]);
					break;
				}
			}
		}
		
		// Calcolo dell'area della faccia
		int N = orderedVertices.size();
		if (N < 3) {
			cerr << "Warning: faccia duale con meno di 3 vertici (id: " << v << ")" << endl;
		} 
		else {
			// Calcolo baricentro della faccia
			Vector3d C = Vector3d::Zero();
			for (int idv : orderedVertices){
				C += dual.Cell0DsCoordinates.col(idv);
			}
			C /= N;

			// Calcolo area sommando triangoli (C, vi, vi+1)
			double area = 0.0;
			for (int i = 0; i < N; ++i) {
				Vector3d p0 = dual.Cell0DsCoordinates.col(orderedVertices[i]);
				Vector3d p1 = dual.Cell0DsCoordinates.col(orderedVertices[(i + 1) % N]);
				area += 0.5 * ( (p0 - C).cross(p1 - C) ).norm();
			}
			
			//Controllo che l'area non sia nulla
			if (area < 1e-14) {
				cerr << "Warning: area nulla (" << area << "), faccia duale con id " << v << endl;
				return false; 
			}
		}
		
		// Salvataggio dei dati della faccia duale
		dual.Cell2DsId.push_back(v);
		dual.Cell2DsVertices.push_back(orderedVertices);
		dual.Cell2DsEdges.push_back(orderedEdges);
	}	

		
// Aggiornamento delle celle 3D
dual.Cell3DsVertices.push_back(dual.Cell0DsId);
dual.Cell3DsEdges.push_back(dual.Cell1DsId);
dual.Cell3DsFaces.push_back(dual.Cell2DsId);
return true;

}

}
