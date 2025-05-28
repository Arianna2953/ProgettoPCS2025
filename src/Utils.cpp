#include "Utils.hpp"
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

int CheckAddEdges(PolyhedralMesh& poly, const Vector2i edge, int& id_edge){
	int w0 = edge[0]; 
	int w1 = edge[1]; 
	
	for (unsigned int i = 0; i < poly.Cell1DsId.size(); i++){
		int u0 = poly.Cell1DsExtrema(0,i);
		int u1 = poly.Cell1DsExtrema(1,i);
		
		if((w0 == u0 && w1 == u1)||(w0 == u1 && w1 == u0)){
			return i;//id del lato che verrebbe duplicato
		}
	}
	id_edge++;
	poly.Cell1DsId.push_back(id_edge);
	poly.Cell1DsExtrema.col(id_edge) = edge;
	return id_edge;	
}

void TriangulationTypeI(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& p, const int& q,const int& n){
	
	using namespace PolyhedralLibrary;
	
	int idV_new = 0;// id dei nuovi vertici del poliedro generato dopo la triangolazione
	int idE_new = -1;// id dei nuovi lati del poliedro generato dopo la triangolazione
	int idF_new = 0;// id delle nuove facce del poliedro generato dopo la triangolazione
	
	int T = n*n;
	int V,E,F;
	
	//relazioni fornite riguardo il numero di vertici, lati e facce del poliedro geodetico generato (presi p, q, b ,c in input
		
	if(p == 3 && q == 3){
		V = 2*T+2;
		E = 6*T;
		F = 4*T;
	}
	else if(p == 3 && q == 4){
		V = 4*T +2;
		E = 12*T;
		F = 8*T;
	}
	else if(p == 3 && q == 5){
		V = 10*T+2;
		E = 30*T;
		F = 20*T;
	}
	
	//predispongo la struttura della PolyhedralMesh andando ad allocare sufficiente memoria per tutte le strutture dati utilizzate per la memorizzazione delle informazioni
	polyNew.NumCell0Ds = V;
	polyNew.Cell0DsId.reserve(V);
	polyNew.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, V);
	
	polyNew.NumCell1Ds = E;
	polyNew.Cell1DsId.reserve(E);
	polyNew.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, E);
	
	polyNew.NumCell2Ds = F;
	polyNew.Cell2DsId.reserve(F);
	polyNew.Cell2DsEdges.reserve(F);
	polyNew.Cell2DsVertices.reserve(F);
	
	polyNew.NumCell3Ds = 1;
	polyNew.Cell3DsId.reserve(1);
	polyNew.Cell3DsEdges.reserve(1);
	polyNew.Cell3DsVertices.reserve(1);
	polyNew.Cell3DsFaces.reserve(1);
		
	for (int f = 0; f < polyOld.NumCell2Ds; f++) {//itero sulle facce del poliedro di base
		
		//mi salvo le coordinate dei vertici della faccia del poliedro di partenza che voglio triangolare
		Vector3d v0 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][0]);
		Vector3d v1 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][1]);
		Vector3d v2 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][2]);
		
		//genero i nuovi vertici, assegnando a ognuno una diversa combinazione di coefficienti (coordinate baricentriche)
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		Vector3d new_point(0.0, 0.0, 0.0);
		MatrixXi vertPosition = MatrixXi::Zero(n+1, n+1);
		for(int i=0;i<=n;i++){ //itero per creare tutte le combinazioni di coefficienti possibili 
            for (int j = 0; j<=i;j++){				
                a=double(i-j)/n;
                b=double(j)/n;
                c=1.0-a-b;
					
				new_point = a*v0 + b*v1 + c*v2;
				//proietto il nuovo vertice creato sulla superficie di una sfera unitaria centrata nell'origine
				if (new_point.norm() < 1e-16) {
					cerr << "Warning: il vettore considerato ha lunghezza nulla";
					break;}
				//new_point.normalize();
				
				//se un verice è già presente, non lo reinserisco
				bool is_copied = false;		
				for(int k = 0; k < idV_new;k++){
					if((polyNew.Cell0DsCoordinates.col(k) - new_point).norm() < 1e-12){
						is_copied = true;
						vertPosition(i,j) = k;//dovrebbe inserire al posto del nuovo id, quello del corrispondente già creato
						break;
					}
				}
				if(!is_copied){
					vertPosition(i,j) = idV_new;//inserisco nella matrice per la "posizione"
					polyNew.Cell0DsId.push_back(idV_new);
					
					
					polyNew.Cell0DsCoordinates.col(idV_new) = new_point;
					
					idV_new++;
				}
			}
		}
		//genero i triangoli (nuovi lati e facce)
		for(int i=0;i < n;i++){
				for (int j = 0; j<=i;j++){
					int vert1 = vertPosition(i,j);
					int vert2 = vertPosition(i+1,j);
					int vert3 = vertPosition(i+1,j+1);
					polyNew.Cell2DsVertices.push_back({vert1,vert2,vert3});
					
					int id1 = CheckAddEdges(polyNew,{vert1,vert2},idE_new);//{vert1,vert2} corrisponde al primo lato
					int id2 = CheckAddEdges(polyNew,{vert2,vert3},idE_new);
					int id3 = CheckAddEdges(polyNew,{vert3,vert1},idE_new);				
					
					polyNew.Cell2DsEdges.push_back({id1,id2,id3});
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
					
					if(i!=j){
						int vert4 = vertPosition(i,j+1);
						polyNew.Cell2DsVertices.push_back({vert1,vert3,vert4});
						
						int id1 = CheckAddEdges(polyNew,{vert1,vert3},idE_new);			
						int id2 = CheckAddEdges(polyNew,{vert3,vert4},idE_new);
						int id3 = CheckAddEdges(polyNew,{vert4,vert1},idE_new);						
						
						polyNew.Cell2DsEdges.push_back({id1,id2,id3});
					
						polyNew.Cell2DsId.push_back(idF_new);
						idF_new++;
					}					
				}
		}		
	}
	
	polyNew.Cell3DsId.push_back(0); // ho solo un poliedro, con identificativo 0 (?)
	polyNew.Cell3DsVertices.push_back(polyNew.Cell0DsId);
	polyNew.Cell3DsEdges.push_back(polyNew.Cell1DsId);
	polyNew.Cell3DsFaces.push_back(polyNew.Cell2DsId);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

void TriangulateFaces(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew,const int& p,const int& q,const int& b,const int& c)
{
	//ho impostato solo un controllo sui valori, bisogna implementare ancora tutta la funzione per la triangolazione!!! 
	//imposto i vari casi in base al valore di b e c in input_iterator
	if((b==0 && c >=1) || (b>=1 && c==0)){
		int n = max(b,c);
		//triangolazione tipo 1
		cout << "Triangolazione di 'tipo 1'" << endl;
		TriangulationTypeI(polyOld, polyNew,p,q,n);
	}
	else if(b==c && b!=0){
		int n = b;
		//triangolazione tipo 2
		cout << "Tipo 2" << endl;
	}
	else {
		cout << "valori di b e c non validi" << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

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
}
}
