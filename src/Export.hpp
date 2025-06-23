#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
// Esportazione delle proprietà di Cell0D dalla PolyhedralMesh
// mesh: struttura dati di tipo PolyhedralMesh
// restituisce true se la creazione della nuova mesh ha avuto successo, altrimenti restituisce false
bool Exportfile0Ds(const PolyhedralMesh& polyNew);

// Esportazione delle proprietà di Cell1D dalla PolyhedralMesh
// mesh: struttura dati di tipo PolyhedralMesh
// restituisce true se la creazione della nuova mesh ha avuto successo, altrimenti restituisce false
bool Exportfile1Ds(const PolyhedralMesh& polyNew);

// Esportazione delle proprietà di Cell2D dalla PolyhedralMesh
// mesh: struttura dati di tipo PolyhedralMesh
// restituisce true se la creazione della nuova mesh ha avuto successo, altrimenti restituisce false
bool Exportfile2Ds(const PolyhedralMesh& polyNew);

// Esportazione delle proprietà di Cell3D dalla PolyhedralMesh
// mesh: struttura dati di tipo PolyhedralMesh
// restituisce true se la creazione della nuova mesh ha avuto successo, altrimenti restituisce false
bool Exportfile3Ds(const PolyhedralMesh& polyNew);

// funzione generale per esportare il poliedro dalla PolyhedralMesh
bool ExportPolyhedron(const PolyhedralMesh& polyNew);

}
