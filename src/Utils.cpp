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

bool TriangulateFaces(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew, const int& b,const int& c)
{
	//ho impostato solo un controllo sui valori, bisogna implementare ancora tutta la funzione per la triangolazione!!! 
	
	
	
	//imposto i vari casi in base al valore di b e c in input_iterator
	if((b==0 && c >=1) || (b>=1 && c==0)){
		int n = max(b,c);
		//triangolazione tipo 1
		cout << "Tipo 1" << endl;
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