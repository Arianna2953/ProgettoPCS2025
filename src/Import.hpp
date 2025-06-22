#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{
// Importa la polyhedral mesh e ne verifica la correttezza
// mesh: struttura dati di tipo PolyhedralMesh da riempire
// file0Ds: percorso del file contenente le informazioni sulle Cell0Ds
// file1Ds: percorso del file contenente le informazioni sulle Cell1Ds
// file2Ds: percorso del file contenente le informazioni sulle Cell2Ds
// file3Ds: percorso del file contenente le informazioni sulle Cell3Ds
// restituisce true se l'importazione ha avuto successo, false altrimenti
bool ImportMesh(PolyhedralMesh& mesh,const string& file0Ds,const string& file1Ds,const string& file2Ds,const string& file3Ds);

// Importa le proprietà delle Cell0Ds dal file Cell0Ds.txt
// file0Ds: percorso del file contenente le informazioni sulle Cell0Ds
// mesh: struttura dati di tipo PolyhedralMesh da aggiornare
// Restituisce true se l'importazione ha avuto successo, false altrimenti
bool ImportCell0Ds(const string& file0Ds, PolyhedralMesh& mesh);

// Importa le proprietà delle Cell0Ds dal file Cell1Ds.txt
// file1Ds: percorso del file contenente le informazioni sulle Cell1Ds
// mesh: struttura dati di tipo PolyhedralMesh da aggiornare
// Restituisce true se l'importazione ha avuto successo, false altrimenti
bool ImportCell1Ds(const string& file1Ds, PolyhedralMesh& mesh);

// Importa le proprietà delle Cell0Ds dal file Cell2Ds.txt
// file2Ds: percorso del file contenente le informazioni sulle Cell2Ds
// mesh: struttura dati di tipo PolyhedralMesh da aggiornare
// Restituisce true se l'importazione ha avuto successo, false altrimenti
bool ImportCell2Ds(const string& file2Ds, PolyhedralMesh& mesh);

// Importa le proprietà delle Cell3Ds dal file Cell3Ds.txt
// file3Ds: percorso del file contenente le informazioni sulle Cell3Ds
// mesh: struttura dati di tipo PolyhedralMesh da aggiornare
// Restituisce true se l'importazione ha avuto successo, false altrimenti
bool ImportCell3Ds(const string& file3Ds, PolyhedralMesh& mesh);

}
