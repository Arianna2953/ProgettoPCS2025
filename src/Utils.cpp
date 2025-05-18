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
void TriangulationTypeI(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& n){
	
	using namespace PolyhedralLibrary;
	
	unsigned int idV_new = 0;// id dei vertici del poliedro generato dopo la triangolazione
	unsigned int idE_new = 0;// id dei lati del poliedro generato dopo la triangolazione
	unsigned int idF_new = 0;// id delle facce del poliedro generato dopo la triangolazione
	unsigned int idV_copy = 0;// id dei vertici del poliedro duplicati (su spigoli o vertici)
	unsigned int idE_copy = 0;// id dei lati del poliedro (sugli spigoli o tra i triangoli della stessa faccia)
	
	unsigned int tot_vertices = (polyOld.NumCell2Ds)*((n+1)*(n+2)/2);
	unsigned int tot_edges = 3*n*n*(polyOld.NumCell2Ds);
	unsigned int tot_faces = n*n*(polyOld.NumCell2Ds);
	
	polyNew.NumCell0Ds = tot_vertices;
	polyNew.Cell0DsId.reserve(tot_vertices);
	polyNew.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, tot_vertices);
	
	polyNew.NumCell1Ds = tot_edges;
	polyNew.Cell1DsId.reserve(tot_edges);
	polyNew.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, tot_edges);
	
	polyNew.NumCell2Ds = tot_faces;
	polyNew.Cell2DsId.reserve(tot_faces);
	polyNew.Cell2DsEdges.reserve(tot_faces);
	polyNew.Cell2DsVertices.reserve(tot_faces);
		
	for (unsigned int f = 0; f < polyOld.NumCell2Ds; f++) {//itero sulle facce del poliedro di base
		
		//cout << "faccia = " << f << endl;
		
		Vector3d v0 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][0]);
		Vector3d v1 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][1]);
		Vector3d v2 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][2]);

		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		Vector3d new_point(0.0, 0.0, 0.0);
		MatrixXd vertPosition = MatrixXd::Zero(n+1, n+1);
		for(unsigned int i=0;i<=n;i++){ //itero per creare tutte le combinazioni di coefficienti possibili 
            for (unsigned int j = 0; j<=i;j++){
				unsigned int s = i-j;
					
                a=double(s)/n;
                b=double(j)/n;
                c=1.0-a-b;
					
				new_point = a*v0 + b*v1 + c*v2;
				
				//se un verice è già presente, non lo reinserisco
				bool is_copied = false;		
				for(int k = 0; k < idV_new;k++){
					if((polyNew.Cell0DsCoordinates.col(k) - new_point).norm() < 1e-16){
						is_copied = true;
						idV_copy = k;//dovrebbe inserire al posto del nuovo id, quello del corrispondente già creato
						break;
					}
				}
				if(is_copied){
					vertPosition(i,j) = idV_copy;
				}
				else {
					if (polyNew.Cell0DsCoordinates.cols() <= idV_new) {
						polyNew.Cell0DsCoordinates.conservativeResize(3, idV_new + 1);
					}
					vertPosition(i,j) = idV_new;//inserisco nella matrice per la "posizione"
					polyNew.Cell0DsId.push_back(idV_new);
					//projectOntoUnitSphere(new_point);
					polyNew.Cell0DsCoordinates.col(idV_new) = new_point;
					
					idV_new++;
				}
			}
		}
		//genero i triangoli (nuovi lati e facce)
		for(unsigned int i=0;i < n;i++){
				for (unsigned int j = 0; j<=i;j++){
					polyNew.Cell1DsExtrema(0,idE_new) = vertPosition(i,j);//v0:polyNew.Cell1DsExtrema(0,idE_new)
					polyNew.Cell1DsExtrema(1,idE_new) = vertPosition(i+1,j);//v1:polyNew.Cell1DsExtrema(1,idE_new)
					idE_new++;
					polyNew.Cell1DsExtrema(0,idE_new) = vertPosition(i+1,j);
					polyNew.Cell1DsExtrema(1,idE_new) = vertPosition(i+1,j+1);
					idE_new++;
					polyNew.Cell1DsExtrema(0,idE_new) = vertPosition(i+1,j+1);
					polyNew.Cell1DsExtrema(1,idE_new) = vertPosition(i,j);
					idE_new++;
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
					
					if(i!=j){
						polyNew.Cell1DsExtrema(0,idE_new) = vertPosition(i,j);
						polyNew.Cell1DsExtrema(1,idE_new) = vertPosition(i+1,j+1);
						idE_new++;
						polyNew.Cell1DsExtrema(0,idE_new) = vertPosition(i+1,j+1);
						polyNew.Cell1DsExtrema(1,idE_new) = vertPosition(i,j+1);
						idE_new++;
						polyNew.Cell1DsExtrema(0,idE_new) = vertPosition(i,j+1);
						polyNew.Cell1DsExtrema(1,idE_new) = vertPosition(i,j);
						idE_new++;
						
						polyNew.Cell2DsId.push_back(idF_new);
						idF_new++;
					}					
				}
		}
		
		
		
		
	}
	//cout << "Dimensione matrice: "
    //      << polyNew.Cell0DsCoordinates.rows() << " x " << polyNew.Cell0DsCoordinates.cols() << std::endl;
		  
	polyNew.Cell0DsCoordinates.conservativeResize(Eigen::NoChange, idV_new);
	//cout << "Dimensione nuova matrice ridotta (no duplicati): "
    //      << polyNew.Cell0DsCoordinates.rows() << " x " << polyNew.Cell0DsCoordinates.cols() << std::endl;
	

}


	
void TriangulateFaces(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew,const int& b,const int& c)
{
	//ho impostato solo un controllo sui valori, bisogna implementare ancora tutta la funzione per la triangolazione!!! 
	
	
	
	//imposto i vari casi in base al valore di b e c in input_iterator
	if((b==0 && c >=1) || (b>=1 && c==0)){
		int n = max(b,c);
		//triangolazione tipo 1
		cout << "Tipo 1" << endl;
		TriangulationTypeI(polyOld, polyNew,n);
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

bool projectOntoUnitSphere(Vector3d& v){
	double len = v.norm();//sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z())
    if (len < 1e-16) {
        cerr << "Warning: il vettore considerato ha lunghezza nulla";
        return false;
    }
    v /= len;
	return true;
}

}