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
	for (unsigned int i = 0; i < poly.Cell1DsId.size(); i++){
		int u0 = poly.Cell1DsExtrema(0,i);
		int u1 = poly.Cell1DsExtrema(1,i);
		int w0 = edge[0]; 
		int w1 = edge[1]; 
		
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
					
					//proietto il nuovo vertice creato sulla superficie di una sfera unitaria centrata nell'origine
					if (new_point.norm() < 1e-16) {
						cerr << "Warning: il vettore considerato ha lunghezza nulla";
						break;}
					new_point.normalize();
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
}
