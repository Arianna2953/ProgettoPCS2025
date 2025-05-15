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
void TriangulationTypeI(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& p,const int& q, const int& n){
	
	using namespace PolyhedralLibrary;
	
	unsigned int T = n*n;
	unsigned int V,E,F;
	
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
	unsigned int id_new = 0;
	polyNew.NumCell0Ds = V; //numero dei nuovi vertici del poliedro
	polyNew.Cell0DsId.reserve(V);
	polyNew.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, 2000);
	//da sistemare le dimensioni!!!!!!
	
	
	polyNew.NumCell1Ds = E;//numero dei nuovi lati del poliedro
	polyNew.Cell1DsId.reserve(E);
	polyNew.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, E);
	
	polyNew.Cell2DsId.reserve(F);
    polyNew.Cell2DsVertices.reserve(F);
    polyNew.Cell2DsEdges.reserve(F);
	
	unsigned int new_id = 0;
	
	for (unsigned int f = 0; f < polyOld.NumCell2Ds; f++) {
		unsigned int faceId = polyOld.Cell2DsId[f];
		const vector<unsigned int>& faceVertices = polyOld.Cell2DsVertices[f];
		const vector<unsigned int>& faceEdges = polyOld.Cell2DsEdges[f];


		//per visualizzazione delle info sulle facce (da togliere poi)
		/*cout << "Faccia ID: " << faceId << endl;

		cout << "  Vertici (ID): ";
		for (unsigned int vId : faceVertices) {
			cout << vId << " ";
		}
		cout << endl;

		cout << "  Coordinate vertici:" << endl;
		for (unsigned int vId : faceVertices) {
			if (vId >= polyOld.Cell0DsCoordinates.cols()) {
				cerr << "    [ERRORE: ID vertice " << vId << " fuori dal range!]" << endl;
				continue;
			}

			Eigen::Vector3d coord = polyOld.Cell0DsCoordinates.col(vId);
			cout << "    ID " << vId << ": " << coord.transpose() << endl;
		}

		cout << "  Lati (ID): ";
		for (unsigned int eId : faceEdges) {
			cout << eId << " ";
		}
		cout << endl << "------------------------------" << endl;*/
		
		
		// Recupero coordinate del primo vertice
		unsigned int id0 = faceVertices[0];
		auto it0 = find(polyOld.Cell0DsId.begin(), polyOld.Cell0DsId.end(), id0);
		if (it0 == polyOld.Cell0DsId.end()) continue;
		unsigned int idx0 = distance(polyOld.Cell0DsId.begin(), it0);
		double x0 = polyOld.Cell0DsCoordinates(0, idx0);
		double y0 = polyOld.Cell0DsCoordinates(1, idx0);
		double z0 = polyOld.Cell0DsCoordinates(2, idx0);
		Vector3d v0(x0, y0, z0);

		// Recupero coordinate del secondo vertice
		unsigned int id1 = faceVertices[1];
		auto it1 = find(polyOld.Cell0DsId.begin(), polyOld.Cell0DsId.end(), id1);
		if (it1 == polyOld.Cell0DsId.end()) continue;
		unsigned int idx1 = distance(polyOld.Cell0DsId.begin(), it1);
		double x1 = polyOld.Cell0DsCoordinates(0, idx1);
		double y1 = polyOld.Cell0DsCoordinates(1, idx1);
		double z1 = polyOld.Cell0DsCoordinates(2, idx1);
		Vector3d v1(x1, y1, z1);

		// Recupero coordinate del terzo vertice
		unsigned int id2 = faceVertices[2];
		auto it2 = find(polyOld.Cell0DsId.begin(), polyOld.Cell0DsId.end(), id2);
		if (it2 == polyOld.Cell0DsId.end()) continue;
		unsigned int idx2 = distance(polyOld.Cell0DsId.begin(), it2);
		double x2 = polyOld.Cell0DsCoordinates(0, idx2);
		double y2 = polyOld.Cell0DsCoordinates(1, idx2);
		double z2 = polyOld.Cell0DsCoordinates(2, idx2);
		Vector3d v2(x2, y2, z2);

		// Ora ho v0, v1, v2 disponibili
		//cout << "Faccia " << polyOld.Cell2DsId[f] << ":\n";
		/*cout << "  v0 = " << v0.transpose() << "\n";
		cout << "  v1 = " << v1.transpose() << "\n";
		cout << "  v2 = " << v2.transpose() << "\n";
		*/
		
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		Eigen::Vector3d new_point(0.0, 0.0, 0.0);
		for(unsigned int i=0;i<=n;i++){ //itero per creare tutte le combinazioni di coefficienti possibili 
              for (unsigned int j = 0; j<=i;j++){
					unsigned int s = i-j;
					
                    a=static_cast<double>(s)/n;
                    b=static_cast<double>(j)/n;
                    c=1.0-a-b;
					
					int marker = 0;  // valore di default

					int zero_count = 0;
					if (a == 0.0) zero_count++;
					if (b == 0.0) zero_count++;
					if (c == 0.0) zero_count++;

					if (zero_count == 2) {
						marker = 2;  // due zeri: vertice
					} else if (zero_count == 1) {
						marker = 1;  // un solo zero: spigolo
					}
				
					//cout << "a=" << a << ",b=" << b << ", c= " << c<<endl;
					new_point = a*v0 + b*v1 + c*v2;
					cout << "  newpoint = " << new_point.transpose() << "\n";
					//controllare
					
					polyNew.Cell0DsId.push_back(id_new);
					polyNew.Cell0DsCoordinates.col(id_new) = new_point;
				
					id_new++;
					//c'è il punto nell'origine... come ci è finito?! :(
					//non funziona con icosaedro (controllo le dimensioni...)
					//bisogna ancora inserire il marker nella struttura, ma ho già messo il controllo
					
			  }
		}
		
	}

}
	

bool TriangulateFaces(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& p,const int& q,const int& b,const int& c)
{
	//ho impostato solo un controllo sui valori, bisogna implementare ancora tutta la funzione per la triangolazione!!! 
	
	
	
	//imposto i vari casi in base al valore di b e c in input_iterator
	if((b==0 && c >=1) || (b>=1 && c==0)){
		int n = max(b,c);
		//triangolazione tipo 1
		cout << "Tipo 1" << endl;
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