
#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
//Calcola il poliedro duale a partire dal poliedro di partenza in 3 step:
// 1. trova i vertici del duale attraverso i baricentri delle facce del poliedro di partenza;
// 2. trova i lati del duale "connettendo" i baricentri di facce adiacenti nel poliedro di partenza;
// 3. trova le facce del duale a partire dai vertici del poliedro di partenza.
//polyhedron: PolyhedralMesh struct in cui sono salvati i dati del poliedro di partenza.
//dual: PolyhedralMesh struct in cui vengono salvati i dati del poliedro duale.
bool DualConstructor(const PolyhedralMesh& polyhedron, PolyhedralMesh& dual);

}
