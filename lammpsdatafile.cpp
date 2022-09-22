#include "lammpsdatafile.h"
#include "molecularmanipulator.h"
#include "utility.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

/* Can only understand LAMMPS data file in the "full" atom style
 * * In every LAMMPS data file, the following info can be read:
* A comment line (including timestep) eg. LAMMPS data file via write_data,
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
    ifstream ifs(filename);
    if(!ifs)
        ERROR("Can't open ["+filename+"] to read.");
    vector<string> lines;
    string line;
    int lineno=0;

    int nAtoms;
    map<int,string> atomElementOfEachType; //inferred from the mass of each atom type

    while(getline(ifs,line)){
        ++lineno;
        lines.push_back(line);
        if(lineno>MY_FILE_LINES_UPPER_BOUND)
            ERROR("File ["+filename+"] too large (>"+
            to_string(MY_FILE_LINES_UPPER_BOUND)+") lines, reading aborted!");
    }
    try {
        // For format of LAMMPS file, especially the header, are not always regular.
        // But the info we look for should always appear within the first 50 lines.
        if(JumpToLine(lines,"atoms",lineno,0,50))
            nAtoms = stoi(StringSplit(lines[lineno])[0]);

        //Read boundary
        JumpToLine(lines,"xlo xhi",lineno,0,lines.size());
        double lohi[3][2];
        for (int i = 0; i < 3; i++) {
            auto parts = StringSplit(lines[lineno + i]);
            for (int j = 0; j < 2; j++) {
                lohi[i][j] = stof(parts[j]);
            }
        }
        JumpToLine(lines,"Masses",lineno,0,lines.size());
        lineno += 2;
        while (StringRegexMatch(lines[lineno], "[0-9]+ [.0-9]+")) {
            auto parts = StringSplit(lines[lineno++]);
            int type = stoi(parts[0]);
            double mass = stof(parts[1]);
            string element = PeriodicTable::PossibleElementWithGivenWeight(mass);
            atomElementOfEachType[type] = element;
        }
//        for (auto &item: atomElementOfEachType)
//            cout << item.first << " " << item.second << endl;
        // Pari Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs, Improper Coeffs are skipped
        // Read Atoms
        JumpToLine(lines,"Atoms",lineno,0,lines.size());
        lineno+=2;
        // The atoms will come in random order. We'll add all atoms into one molecule, then sort_atoms the atoms, and
        // separate the atoms into multiple molecules;
        ms.AddMolecule(Molecule());
        map<string,string> atomSerialToMolSerialMap;
        int nTotalAtomsCount = 0;
        while(StringRegexMatch(lines[lineno],"[0-9]+ [0-9]+ ")){
            // One line for an atom: (the last three numbers are images flags and are ignored)
            // 5 1 8 2.1 2.94997 2.96969 4.25362 0 0 0
            Atom a;
            string molSerial;
            istringstream iss(lines[lineno]);
            iss>>a.globalSerial>>molSerial>>a.type>>a.charge>>a.xyz[0]>>a.xyz[1]>>a.xyz[2];
            a.element = a.name = atomElementOfEachType[stoi(a.type)];
            ms[0].AddAtom(a);
            atomSerialToMolSerialMap[a.globalSerial] = molSerial;
            nTotalAtomsCount++;
            ++lineno;
        }
        // Now sort_atoms the atoms. Note the lambda function compares the globalSerials (converts to int) of two atoms.
        sort(ms[0].atoms.begin(), ms[0].atoms.end(),
             [](auto &a1, auto &a2) {
                    return stoi(a1->globalSerial)<stoi(a2->globalSerial);}
             );
        if(nTotalAtomsCount != nAtoms) { // if not true, the declared # of atoms != actual read atoms.
            ostringstream oss;
            oss << "\nAcutally read atoms (" << nTotalAtomsCount << ") does match declared atom count (" << nAtoms
                << ")";
            output(oss.str());
            throw exception();
        }
        vector<int> map_for_reconstructing(nTotalAtomsCount); // needed for MolSysReorganize() to redistribute atoms into multiple mols;
        for(auto &item:atomSerialToMolSerialMap){
            int atomIndex = stoi(item.first)-1;
            int molIndex = stoi(item.second)-1;
            map_for_reconstructing[atomIndex] = molIndex;
        }
        MolSysReorganize(ms,map_for_reconstructing);

        // "Velocities" are ignored
        // Now we read the bonds. Before that, let's renumber atoms ans make a MolecularAccessor
        ms.RenumberAtomSerials();
        MolecularSystemAccessor msa(ms);
        JumpToLine(lines,"Bonds",lineno,0,lines.size());
        lineno+=2;
        while(StringRegexMatch(lines[lineno],"[0-9]+[\\s]+[0-9]+")){
            Bond b;
            string tmp;
            istringstream iss(lines[lineno]);
            iss>>tmp>>b.type>>b.atom1>>b.atom2;
            auto a1_id = msa.MolAndLocalIndexOfAtom(b.atom1);
            auto a2_id = msa.MolAndLocalIndexOfAtom(b.atom2);
            if(a1_id.first == a2_id.first) {
                // In-molecule bond
                b.atom1 = ms[a1_id.first][a1_id.second].serial;
                b.atom2 = ms[a2_id.first][a2_id.second].serial;
                ms[a1_id.first].bonds.push_back(make_shared<Bond>(b));
            }else{
                // Inter-molecular
                ms.interMolecularBonds.push_back(make_shared<Bond>(b));
            }
            lineno++;
        }
        // No renumbering needed for the second time.
        // Angles, Dihedrals, Impropers are ignored.
        return nTotalAtomsCount>0;
    }catch(exception e){
        ERROR("\nWhile reading line "+to_string(lineno+1)+" of ["+filename+"]:\n\""+line+"\"");
    }
    return false; // should not reach here if no exception
}
bool LAMMPSDataFile::Write(MolecularSystem &ms, string filename){
    // Maybe needs a RTTI to determine whether can be written as LAMMPSData
    ERROR("Sorry, Writing as LAMMPSData file not supported.");
    return false;
};
