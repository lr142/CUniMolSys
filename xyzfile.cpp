#include "xyzfile.h"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

bool XYZFile::Read(MolecularSystem &ms, std::string filename) {
    /* # XYZ file usually contains only element and x,y,and z coordinates.
    # In our system we shall extend the XYZ format a little bit, to allow a trailing flag "T" or "F" to
    # denote whether an atom is flexible in the QM or MD simulation.
    # Note that an XYZ file can contain more than 1 molecule. It's also optional to supply a molecule name
    # in the comment line ( the 2nd line, or the line below each atom count if multiple molecules are present) */
    ifstream ifs(filename);
    if(!ifs)
        ERROR("Can't open xyz file["+filename+"] to read!");
    string line;

    ms.Clear();

    while(getline(ifs,line)){
        int atomCount;
        try {
            atomCount = stoi(StringRemoveComment(line));
        }catch(std::exception& e){
            continue;
        }
        shared_ptr<Molecule> currentMol = make_shared<Molecule>();
        getline(ifs,line);
        currentMol->name = StringStrip(line);
        if(ms.GetName()=="")
            ms.SetName(currentMol->name);
        for(int i=0;i<atomCount;i++) {
            getline(ifs, line);
            auto parts = StringSplit(line);
            if (parts.size() < 4)
                ERROR("Format error in xyz file[" + filename + "]:\n" + line);
            Atom a;
            a.element = parts[0];
            a.xyz = {stof(parts[1]), stof(parts[2]), stof(parts[3])};
            if (parts.size() > 4)
                a.flexible = StringToUpper(parts[4]) == "T" ? true : false;
            currentMol->AddAtom(a);
        }
        ms.molecules.push_back(currentMol);
    }
    ms.RenumberAtomSerials();
    return true;
}

bool XYZFile::Write(MolecularSystem &ms, std::string filename) {
    ofstream ofs(filename);
    if(!ofs)
        ERROR("Can't open xyz file["+filename+"] to write!");
    for(int i=0;i<ms.MoleculesCount();i++){
        ofs<<ms[i].AtomsCount()<<endl;
        ofs<<(ms[i].name==""?ms.GetName():ms[i].name)<<endl;
        for(int j=0;j<ms[i].AtomsCount();j++){
            ofs<<ms[i][j].ShowAsXYZ()<<endl;
        }
    }
    return true;
}
