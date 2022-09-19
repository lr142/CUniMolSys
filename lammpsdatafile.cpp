#include "lammpsdatafile.h"
#include <iostream>
#include <fstream>
using namespace std;

/* Can only understand LAMMPS data file in the "full" atom style
 * * In every LAMMPS data file, the following info can be read:
* A comment line (including timestep) eg. LAMMPS data file via write_data, version 3 Mar 2020, timestep = 10000
* number of atoms, bonds, angle, dihedrals, etc.
* boundary info, e.g:
* -2.0000000000000000e+01 1.7000000000000000e+02 ylo yhi
* -2.0000000000000000e+01 1.7000000000000000e+02 zlo zhi
* Then the following sections:
* Masses
* Pair Coeffs * lj/cut/coul/long
* Bond Coeffs # harmonic
* Angle Coeffs # harmonic
* Dihedral Coeffs # opls
* Improper Coeffs # harmonic
* Atoms # full
* Velocities
* Bonds
* Angles
* Not all sections may be read. For example velocities may be absent.
* In this function we shall read in only boundary, masses, atoms, bonds.


* No element type is provided, but can be inferred from masses.
* Some atomic weights may not be included in the list,
* such as some united-atom. The user need to modify them manually later.
* An massless (<0.1) atom is recognized as pseudo atom 'M'.
 */

bool LAMMPSDataFile::Read(MolecularSystem &ms, string filename){
#define REPORT_ERR ERROR("\nWhile reading line "+to_string(lineno+1)+" of ["+filename+"]:\n\""+line+"\"")

    ifstream ifs(filename);
    if(!ifs)
        ERROR("Can't open ["+filename+"] to read.");
    vector<string> lines;
    string line;
    int lineno=0;

    int nAtoms;
    int nAtomTypes;
    int nBonds;
    int nBondTypes;

    while(getline(ifs,line)) {
        ++lineno;
        if (!StringEndsWith(line, "atom", false))
            continue;
        else
            nAtoms = stoi(StringSplit(line)[0]);
    }
    getline(ifs,line);
    nAtomTypes = stoi(StringSplit(line)[0]);
    getline(ifs,line);
    nBonds = stoi(StringSplit(line)[0]);
    getline(ifs,line);
    nBondTypes = stoi(StringSplit(line)[0]);
    lineno+= 3;




    }
    return true;
}
bool LAMMPSDataFile::Write(MolecularSystem &ms, string filename){
    // Maybe needs a RTTI to determine whether can be written as LAMMPSData
    ERROR("Sorry, Writing as LAMMPSData file not supported.");
    return false;
};
