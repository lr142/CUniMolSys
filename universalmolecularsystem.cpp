#include "universalmolecularsystem.h"
#include <sstream>
using namespace std;

string Atom::ShowAsXYZ(){
    ostringstream oss;
    oss<<element<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2];
    return oss.str();
}
string Atom::ShowAllFields(){
    ostringstream oss;
    oss<<serial<<" "<<element<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2] \
    <<name<<" "<<type<<" "<<flexible<<" "<<charge;
    return oss.str();
}
// returns ith atom. Supports negative index(counting from end)
Atom& Molecule::operator[](int index) {
    int nAtoms = AtomsCount();
    if(index < -nAtoms || index>=nAtoms)
        throw out_of_range("Molecules::[]");
    return *atoms[(index+nAtoms)%nAtoms];
}
Bond& Molecule::GetBond(int index){
    int nBonds = BondsCount();
    if(index<-nBonds || index>=nBonds)
        throw out_of_range("Molecules::getBond");
    return *bonds[(index+nBonds)%nBonds];
}
vector<vector<int>>& Molecule::GetBondedMap(bool update){
    static vector<vector<int>> bondedTo;
    if(!update)
        return bondedTo;
    int nAtoms = AtomsCount();
    bondedTo.resize(nAtoms);
    map<string,int> serialToIndexMap;
    for(int i=0;i<nAtoms;i++) {
        bondedTo[i].clear();
        serialToIndexMap[atoms[i]->serial] = i;
    }
    for(int i=0;i<bonds.size();i++){
        string from = bonds[i]->atom1;
        string to = bonds[i]->atom2;
        auto not_found_pos = serialToIndexMap.end();
        if(serialToIndexMap.find(from)==not_found_pos || serialToIndexMap.find(to)==not_found_pos){
            error("Dangling bond detected: ["+from+"-"+to+"] in molecule serial ="+ serial,false);
        }
        int fromIndex = serialToIndexMap[from];
        int toIndex   = serialToIndexMap[to];
        bondedTo[fromIndex].push_back(toIndex);
        bondedTo[toIndex].push_back(fromIndex);
    }
    return bondedTo;
}
string Molecule::Summary() {
    ostringstream oss;
    oss<<"Molecule Serial: "<<serial<<", Name: "<<name<<", Type: "<<type
    <<", with "<<AtomsCount()<<" Atoms, and "<<BondsCount()<<" Bonds";
    return oss.str();
}

Molecule Molecule::DeepCopy(){
    // Firstly, call the operator =, which returns a shallow copy.
    Molecule copy = *this;
    // Then need to copy atoms and bonds:
    copy.atoms.clear();
    for(int i=0;i<this->AtomsCount();i++){
        auto pNewAtom = make_shared<Atom>();
        *pNewAtom = *(this->atoms[i]);
        copy.atoms.push_back(pNewAtom);
    }
    copy.bonds.clear();
    for(int i=0;i<this->BondsCount();i++){
        auto pNewBond = make_shared<Bond>();
        *pNewBond = *(this->bonds[i]);
        copy.bonds.push_back(pNewBond);
    }
    return copy;
}

bool Molecule::ConsistencyCheck() {
    // 1. Atom序列号必须唯一
    map<string,int> serialToIndexMap;
    for(int i=0;i<AtomsCount();i++){
        if(serialToIndexMap.find(atoms[i]->serial) != serialToIndexMap.end())
            return false;
        serialToIndexMap[atoms[i]->serial] = i;
    }
    // 2. 没有悬空键
    for(int i=0;i<BondsCount();i++){
        string from = bonds[i]->atom1;
        string to = bonds[i]->atom2;
        if(serialToIndexMap.find(from) == serialToIndexMap.end())
            return false;
        if(serialToIndexMap.find(to) == serialToIndexMap.end())
            return false;
    }
    return true;
}

Boundary::Boundary(bool ortho):orthogonal_(ortho){}
void Boundary::SetUVW(XYZ uvw[3]){
    this->uvw[0] = uvw[0];
    this->uvw[1] = uvw[1];
    this->uvw[2] = uvw[2];
}
void Boundary::SetUVW(XYZ u,XYZ v,XYZ w){
    uvw[0]=u; uvw[1]=v; uvw[2]=w;
}
void Boundary::SetUVW(double uvw[3][3]) {
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            this->uvw[i][j] = uvw[i][j];
        }
    }
}
void Boundary::SetOrigin(XYZ origin) {
    this->origin = origin;
}
XYZ& Boundary::operator[](int index){
    if(index<0 || index>2)
        throw out_of_range("Boundary::operator[]");
    return this->uvw[index];
}
void Boundary::SetLoHi(double lohi[3][2]){
    if(!Orthogonal())
        throw runtime_error("In Boundary::SetLoHi(), boundary not orthogonal!");
    uvw[0] = {lohi[0][1]-lohi[0][0], 0, 0};
    uvw[1] = {0, lohi[1][1]-lohi[1][0], 0};
    uvw[2] = {0, 0, lohi[2][1]-lohi[2][0]};
    origin = {lohi[0][0], lohi[1][0], lohi[2][0]};
}
void Boundary::GetLoHi(double lohi[3][2]){
    if(!Orthogonal())
        throw runtime_error("In Boundary::GetLoHi(), boundary not orthogonal!");
    for(int i=0;i<3;i++){
        lohi[i][0] = origin[i];
        lohi[i][1] = lohi[i][0] + uvw[i][i];
    }
}
void Boundary::GetUVW(double uvw[3][3]){
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            uvw[i][j] = this->uvw[i][j];
        }
    }
}
bool Boundary::Periodic() {
    return uvw[0]==XYZ(0,0,0) and uvw[1]==XYZ(0,0,0) and uvw[2]==XYZ(0,0,0);
}
string Boundary::Show(){
    ostringstream oss;
    oss<<"UVW: "<<(*this)[0]<<", "<<(*this)[1]<<", "<<(*this)[2]<<endl;
    oss<<"Center: "<<origin<<endl;
    double lohi[3][2];
    GetLoHi(lohi);
    oss<<"LoHi: { ["<<lohi[0][0]<<", "<<lohi[0][1]<<"], ["
    <<lohi[1][0]<<", "<<lohi[1][1]<<"], ["
    <<lohi[2][0]<<", "<<lohi[2][1]<<"] }";
    return oss.str();
}

MolecularSystem::MolecularSystem(string name): name(name){}
void MolecularSystem::Read(MolecularFile *pFile,string filename){
    pFile -> Read(*this,filename); //Strategy Pattern
}
void MolecularSystem::Write(MolecularFile *pFile,string filename){
    pFile -> Write(*this,filename); //Strategy Pattern
}
MolecularSystem MolecularSystem::DeepCopyWithTrajectory(){
    MolecularSystem copy = this->DeepCopy();
    /* not yet implemented
    copy.trajectory = make_shared<Trajectory>()
     */
    return copy;
}
MolecularSystem MolecularSystem::DeepCopy(){
    MolecularSystem copy(this->name);
    for(int i=0;i<MoleculesCount();i++){
        shared_ptr<Molecule> pNewMol = make_shared<Molecule>();
        *pNewMol = this->molecules[i]->DeepCopy();
        copy.molecules.push_back(pNewMol);
    }
    for(int i=0;i<this->interMolecularBonds.size();i++){
        copy.interMolecularBonds.push_back(make_shared<Bond>());
        *(copy.interMolecularBonds[i]) = *(this->interMolecularBonds[i]);
    }
    return copy;
}
int MolecularSystem::AtomsCount(){
    int count = 0;
    for(int i=0;i<MoleculesCount();i++){
        count += (*this)[i].AtomsCount();
    }
    return count;
}
int MolecularSystem::MoleculesCount() {
    return molecules.size();
}
Molecule& MolecularSystem::operator[](int index){
    int molCount = MoleculesCount();
    if(index < -molCount || index >= molCount)
        throw out_of_range("MolecularSystem::operator[]");
    return *molecules[(index+molCount)%molCount];
}
shared_ptr<Atom> MolecularSystem::GetAtom(int index){
    int nAtoms = AtomsCount();
    if(index<-nAtoms || index>=nAtoms)
        return nullptr;
    index = (index+nAtoms)%nAtoms;
    for(int i=0;i<MoleculesCount();i++){
        if(index < molecules[i]->AtomsCount())
            return molecules[i]->atoms[index];
        else
            index -= molecules[i]->AtomsCount();
    }
    return nullptr;
}
int MolecularSystem::SearchMoleculeBySerial(string molSerial,bool update){
    static map<string,int> serialToIndexMap;
    if(serialToIndexMap.empty() || update){
        serialToIndexMap.clear();
        for(int i=0;i<MoleculesCount();i++){
            serialToIndexMap[(*this)[i].serial] = i;
        }
    }
    if(serialToIndexMap.find(molSerial) != serialToIndexMap.end())
        return serialToIndexMap[molSerial];
    else
        return -1;
}
std::pair<int,int> MolecularSystem::SearchAtomByGlobalSerial(string atomGlobalSerial,bool update){
    static map<string,pair<int,int>> serialToIndexMap;
    if(serialToIndexMap.empty() || update){
        serialToIndexMap.clear();
        for(int i=0;i<MoleculesCount();i++){
            for(int j=0;j<(*this)[i].AtomsCount();j++){
                serialToIndexMap[(*this)[i][j].globalSerial] = make_pair(i,j);
            }
        }
    }
    if(serialToIndexMap.find(atomGlobalSerial) != serialToIndexMap.end())
        return serialToIndexMap[atomGlobalSerial];
    else
        return make_pair(-1,-1);
}
void MolecularSystem::RenumberAtomSerials(int startingGlobalSerial) {
    // Call SearchAtomByGlobalSerial(true) to record the old
    // Serials locations.
    SearchAtomByGlobalSerial("",true);
    // Also, need to record the old atoms serials within each molecule
    vector<map<string,int>> oldInMoleculeSerialToIndexMaps;
    for(int iMol=0;iMol<MoleculesCount();iMol++){
        map<string,int> inMolSerialToIndexMap;
        for(int iAtom=0;iAtom<(*this)[iMol].AtomsCount();iAtom++){
            inMolSerialToIndexMap[(*this)[iMol][iAtom].serial] = iAtom;
        }
        oldInMoleculeSerialToIndexMaps.push_back(inMolSerialToIndexMap);
    }
    // Start the renumbering
    int runningMolSerial = startingGlobalSerial;
    int runningAtomGlobalSerial = startingGlobalSerial;
    for(int iMol=0;iMol<MoleculesCount();iMol++){
        (*this)[iMol].serial = to_string(runningMolSerial++);
        int runningAtomLocalSerial = 1;
        for(int iAtom=0;iAtom<(*this)[iMol].AtomsCount();iAtom++){
            (*this)[iMol][iAtom].serial = to_string(runningAtomLocalSerial++);
            (*this)[iMol][iAtom].globalSerial = to_string(runningAtomGlobalSerial++);
        }
        // Now take care of the in-molecule bonds
        for(int iBond=0;iBond<(*this)[iMol].BondsCount();iBond++){
            auto pbond = (*this)[iMol].bonds[iBond];
            int iFrom = oldInMoleculeSerialToIndexMaps[iMol][pbond->atom1];
            pbond->atom1 = (*this)[iMol][iFrom].serial;
            int iTo   = oldInMoleculeSerialToIndexMaps[iMol][pbond->atom2];
            pbond->atom2 = (*this)[iMol][iTo].serial;
        }
    }
    // Finally take care of the inter-molecule bonds:
    for(int iBond=0;iBond<this->interMolecularBonds.size();iBond++){
        auto pBond = this->interMolecularBonds[iBond];
        // Must call the search function by update=false, which gives the old info
        auto iFrom = SearchAtomByGlobalSerial(pBond->atom1,false);
        auto iTo   = SearchAtomByGlobalSerial(pBond->atom2,false);
        auto strFrom = (*this)[iFrom.first][iFrom.second].globalSerial;
        auto strTo     = (*this)[iTo.first][iTo.second].globalSerial;
        pBond->atom1 = strFrom;
        pBond->atom2 = strTo;
    }
    // This last steps are optional, manually update the static cached data within
    // Various functions
    SearchAtomByGlobalSerial("",true);
    SearchMoleculeBySerial("",true);
}
void MolecularSystem::Clear(){
    molecules.clear();
    interMolecularBonds.clear();
    trajectory = nullptr;
}
string MolecularSystem::Summary(){
    int bondCount = 0;
    for(int iMol=0;iMol<MoleculesCount();iMol++){
        bondCount+=(*this)[iMol].BondsCount();
    }
    ostringstream oss;
    oss<<"MolecularSystem "<<name<<" has "<<MoleculesCount()<<" molecules, "
    <<AtomsCount()<<" atoms, "<<bondCount<<" in-molecular bonds, and "
    <<interMolecularBonds.size()<<" inter-molecular bonds";
    return oss.str();
}

void MolecularSystem::Translate(XYZ offset){
    int nAtoms = AtomsCount();
    for(int i=0;i<nAtoms;i++){
        GetAtom(0)->xyz+=offset;
    }
}
void MolecularSystem::Rotate(double clockwise_degree, XYZ axis){
    int nAtom = AtomsCount();
    XYZ* newCoords = new XYZ[nAtom];
    for(int i=0;i<nAtom;i++){
        newCoords[i] = GetAtom(i)->xyz;
    }
    XYZRotate(newCoords, nAtom, clockwise_degree, axis);
    for(int i=0;i<nAtom;i++){
        GetAtom(i)->xyz = newCoords[i] ;
    }
    delete [] newCoords;
}
void MolecularSystem::FractionalToCartesian(){
    if(!boundary.Orthogonal())
        error("MolecularSystem::FractionalToCartesian() not supported for non-orthogonal systems!",true);
    double uvw[3][3];
    boundary.GetUVW(uvw);
    double X,Y,Z;
    X = uvw[0][0];
    Y = uvw[1][1];
    Z = uvw[2][2];
    int nAtom = AtomsCount();
    for(int i=0;i<nAtom;i++){
        auto &xyz = GetAtom(i)->xyz;
        xyz[0] /= X;
        xyz[1] /= Y;
        xyz[2] /= Z;
    }
}
