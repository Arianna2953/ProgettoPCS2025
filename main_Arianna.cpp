#include <iostream>
#include <string>
#include <cctype>  // per isdigit()
#include "PolyhedralMesh.hpp"
#include "ImportExport.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"


using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

PolyhedralMesh regularPolyhedron;
int main(){
Exportfile0Ds(regularPolyhedron);
}