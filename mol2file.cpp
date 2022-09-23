#include <cmath>
#include "mol2file.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
Mol2File::Mol2File(bool writeElementInsteadOfType): MolecularFile(){
    writeElementInsteadOfType_ = writeElementInsteadOfType;
}

bool Mol2File::Read(MolecularSystem &ms, string filename) {
    ifstream ifs(filename);
    if(!ifs)
        ERROR("Can't open ["+filename+"] to read.");
    enum STATE{ATOM,BOND,OTHER,CRYSIN,MOLECULE};
    STATE state = OTHER;

    ms.Clear();
    string line;
    int lineno = 0;
    while(getline(ifs,line)){
        ++lineno;
        if(StringStartsWith(line,"@")){
            // Special Character that means the start of a new section
            if(StringStartsWith(line,"@<TRIPOS>MOLECULE")){
                ms.AddMolecule(Molecule());
                state = MOLECULE;
            }else if(StringStartsWith(line,"@<TRIPOS>ATOM")){
                state = ATOM;
            }else if(StringStartsWith(line,"@<TRIPOS>BOND")){
                state = BOND;
            }else if(StringStartsWith(line,"@<TRIPOS>CRYSIN")) {
                //This is the default way Materials Studio writes PBC boundaries.
                state = CRYSIN;
            }else{
                state = OTHER;
                //more keywords and more features can be read from the MOL2 file. If necessary, extend the program
                //from the above enumeration statements and add respective actions below.
            }
            continue;
        }
        // if reaches here, this line doesn't start with a "@"
        if(state == MOLECULE){
            string molName = StringStrip(line);
            ms[-1].name = molName;
            if(ms.GetName() == "")
                ms.SetName(molName);
            state = OTHER;
        }else if(state==ATOM){
            Atom a;
            if(parseAtomLine(line,a))
                ms[-1].AddAtom(a);
            else
                showError(lineno,line,filename,"",true);
        }else if(state==BOND){
            Bond b;
            if(parseBondLine(line,b))
                ms[-1].bonds.push_back(make_shared<Bond>(b));
            else
                showError(lineno,line,filename,"",true);
        }else if(state==CRYSIN){
            if(!parseBoundaryLine(line, ms.boundary))
                showError(lineno,line,filename,"",true);
        }else{
            // Do nothing
        }
    }

    //Final Check
    if(ms.AtomsCount()==0){
        showError(lineno,"",filename,"No atoms/molecules found in the file.",true);
    }
    ms.RenumberAtomSerials();
    return true;
}


bool Mol2File::parseAtomLine(string line,Atom &a){
    line = StringStrip(line);
    if(line.size()==0)
        return false;
    auto parts = StringSplit(line);
    try{
        /*  # Format:
            # atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]
            # In Protein Structures, the subst_id and subst_name are usually residue serial and residue name, for example:
            # 31 CD2           0.747000   38.794000   38.347000 C.3      4 LEU27      0.0072
            # An additional section called @<TRIPOS>SUBSTRUCTURE is present in the MOL2 file telling that to which residue
            # an atom belongs to. But in this version we ignore those residue information.
            # Just get the charge from parts[8] if it has one.
        */
        a.serial = parts[0];
        a.name = parts[1];
        a.xyz = XYZ(stof(parts[2]),stof(parts[3]),stof(parts[4]));
        a.type = parts[5];
        auto type_splitted = StringSplit(a.type,'.');
        a.element = StringToCapitalized(type_splitted[0]);
        a.charge = parts.size()>8 ? stof(parts[8]) : 0.0;
    }catch(exception){
        return false;
    }
    return true;
}
bool Mol2File::parseBondLine(string line,Bond &b) {
    auto parts = StringSplit(line);
    if(parts.size()<4)
        return false;
    b.atom1 = parts[1];
    b.atom2 = parts[2];
    b.type = parts[3];
    return true;
}
bool Mol2File::parseBoundaryLine(string line, Boundary &boundary){
   // Format: lx ly lz A B C space_group setting
   // Here I'm assuming C is along the Z axis, and A is along the X axis, B is in X-Y plane
   // i.e. A=B=90, C may not to 90.

   auto parts = StringSplit(line);
   XYZ_DTYPE lx,ly,lz,C;
   try{
       lx = stof(parts[0]);
       ly = stof(parts[1]);
       lz = stof(parts[2]);
//       A  = stof(parts[3]);
//       B  = stof(parts[4]);
       C  = stof(parts[5]);
   }catch(exception){
       return false;
   }
   boundary.SetOrigin(XYZ(0,0,0));
   XYZ u = XYZ(lx,0,0);
   XYZ v = XYZ(ly*cos(C/180*MY_PI), ly*sin(C/180*MY_PI),0);
   XYZ w = XYZ(0,0,lz);
   boundary.SetUVW(u,v,w);
   return true;
}
void Mol2File::showError(int lineno,string line,string filename,string comment,bool fatal){
    output("Error while reading ["+filename+"] at line "+ to_string(lineno)+":");
    if(line!="")
        output(line);
    if(comment!="")
        output(comment);
    if(fatal)
        ERROR("Mol2File::Read() failed.");
    else
        WARNING("Mol2File::Read() failed.");
}

bool Mol2File::Write(MolecularSystem &ms, std::string filename) {
    ofstream ofs(filename);
    if(!ofs)
        ERROR("Can't open ["+filename+"] to write.");
    for(int i=0;i<ms.MoleculesCount();i++)
        writeAMolecule(ms[i],ofs);
    if(ms.Periodic()){
        ofs<<"@<TRIPOS>CRYSIN"<<endl;
        XYZ v = ms.boundary.GetV();
        double C = atan(v[1]/v[0])/MY_PI*180.0;
        if(C<0)
            C+=180;
        ofs<<ms.boundary.GetU()[0]<<" "
           <<v.Norm()<<" "
           <<ms.boundary.GetW()[2]<<" "
           <<"90 90 "
           <<C
           <<" 1 1"<<endl;
    }
    return true;
}

void Mol2File::writeAMolecule(Molecule &mol, std::ofstream &ofs) {
    ofs<<"@<TRIPOS>MOLECULE"<<endl;
    ofs<< (mol.name=="" ? "mol" : mol.name)<<endl;
    ofs<<mol.AtomsCount()<<" "<<mol.BondsCount()<<endl;
    ofs<<"SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM"<<endl;
    for(int i=0;i<mol.AtomsCount();i++){
        Atom &a = mol[i];
        string typeSection;
        if( writeElementInsteadOfType_ || a.type=="")
            typeSection = a.element;
        else
            typeSection = a.type;

        ofs<<a.serial<<" ";
        // If you want the name to be unique:
        //ofs<<a.element+a.serial<<" ";
        // Maybe you don't want to rename the atoms. In that case, use the statement below.
        ofs<<a.name<<" ";
        ofs<<a.xyz[0]<<" "<<a.xyz[1]<<" "<<a.xyz[2]<<" "
        <<typeSection<<" "<<"0.0"<<" "<<"****"<<" "<<a.charge<<endl;
    }
    ofs<<"@<TRIPOS>BOND"<<endl;
    for(int i=0;i<mol.BondsCount();i++){
        auto pBond = mol.bonds[i];
        ofs<<i+1<<" "<<pBond->atom1<<" "<<pBond->atom2<<" "<<pBond->type<<endl;
    }
}
