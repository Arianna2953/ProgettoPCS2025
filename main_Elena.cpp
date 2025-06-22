#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"

//#include "ImportExport.hpp"
//#include "Utils.hpp"
#include "Import.hpp"
#include "Export.hpp"
#include "Triangulation.hpp"
#include "Dual.hpp"
#include "ShortPath.hpp"

#include "UCDUtilities.hpp"


using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main() {
    int V = 10;
	vector<double> dist = {10.0, 9.0,  8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
	vector<pair<double,int>> priorityQueue(V);
	for (int i = 0; i < V; i++){
		priorityQueue.push_back(pair(dist[i], i));
	}
	
	cout << "CODA CON PRIORITA" << endl;
	for (auto it = priorityQueue.begin(); it != priorityQueue.end(); it++){
		cout << get<0>(*it) << "  -  " << get<1>(*it) << endl;
	}
	
	//make_heap(priorityQueue.begin(), priorityQueue.end(), greater<pair<double,int>>{});
	make_heap(priorityQueue.begin(), priorityQueue.end());
	
	cout << "CODA CON PRIORITA" << endl;
	for (auto it = priorityQueue.begin(); it != priorityQueue.end(); it++){
		cout << get<0>(*it) << "  -  " << get<1>(*it) << endl;
	}
	
}