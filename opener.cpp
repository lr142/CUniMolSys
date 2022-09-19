#include "opener.h"
#include "mol2file.h"
#include "xyzfile.h"
#include "lammpsdatafile.h"
#include <regex>
using namespace std;



shared_ptr<MolecularFile> FileTypePicker(string filename){
    shared_ptr<MolecularFile> molfile = nullptr;
    if(StringEndsWith(filename,"MOL2",false)){
        molfile = make_shared<Mol2File>();
    }else if(StringEndsWith(filename,"xyz", false)){
        molfile = make_shared<XYZFile>();
    }else if(StringEndsWith(filename,"data",false)){
        molfile = make_shared<LAMMPSDataFile>();
    }else{
        molfile = nullptr;
    }
    return molfile;
}

void QuickOpen(MolecularSystem &ms,string filename){
    shared_ptr<MolecularFile> molfile = FileTypePicker(filename);
    if(molfile!= nullptr){
        ms.Read(molfile.get(),filename);
    }else{
        ERROR("Trying to read :["+filename+"]: File type not understood!");
    }
}

void QuickSave(MolecularSystem &ms,string filename){
    shared_ptr<MolecularFile> molfile = FileTypePicker(filename);
    if(molfile!= nullptr){
        ms.Write(molfile.get(),filename);
    }else{
        ERROR("Trying to write:["+filename+"]: File type not understood!");
    }
}
