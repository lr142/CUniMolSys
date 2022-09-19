#include "opener.h"
#include <regex>
using namespace std;

bool StringEndsWithCaseInsensitive(string str,string pattern){
    pattern = pattern+"$";
    regex r(pattern,regex::icase);
    smatch sresults;
    if(regex_search(str,sresults,r)){
//        cout<<sresults.size()<<endl;
//        cout<<sresults.str()<<endl;
        return true;
    }
    return false;
}
bool StringContainsCaseInsensitive(string str,string pattern){
    regex r(pattern,regex::icase);
    smatch sresults;
    if(regex_search(str,sresults,r)){
        return true;
    }
    return false;
}
shared_ptr<MolecularFile> FileTypePicker(string filename){
    shared_ptr<MolecularFile> molfile = nullptr;
    if(StringEndsWithCaseInsensitive(filename,"MOL2")){
        molfile = make_shared<Mol2File>();
    }else if(StringEndsWithCaseInsensitive(filename,"xyz")){
        molfile = make_shared<XYZFile>();
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
