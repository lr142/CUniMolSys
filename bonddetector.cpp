#include "bonddetector.h"
#include "utility.h"
#include "xyzfile.h"
#include <cmath>
#include <fstream>

using namespace std;

GridForNeighList::GridForNeighList(double SystemLx, double SystemLy, double SystemLz, double gridsize):gridsize_(gridsize) {
    Lx_ = int(ceil(SystemLx/GridSize()));
    Ly_ = int(ceil(SystemLy/GridSize()));
    Lz_ = int(ceil(SystemLz/GridSize()));
    if(Lx_<=0 || Ly_<=0 || Lz_<=0 || gridsize<=0 ){
        error("In GridForNeighList(), one of the dimension <=0");
    }else if(double(Lx_)*Ly_*Lz_ > 1000000000 ){ // > 1 billion, convert to double to prevent overflow
        error("In GridForNeighList(), too many grids: "+to_string(double(Lx_)*Ly_*Lz_)+", check gridsize_");
    }

    grid = new vector<AtomInGrid>** [Lx()];
    for(int i=0;i<Lx();i++){
        grid[i] = new vector<AtomInGrid>* [Ly()];
        for(int j=0;j<Ly();j++){
            grid[i][j] = new vector<AtomInGrid> [Lz()];
        }
    }
}

GridForNeighList::~GridForNeighList(){
    for(int i=0;i<Lx();i++){
        for(int j=0;j<Ly();j++){
            delete [] grid[i][j];
        }
        delete [] grid[i];
    }
    delete [] grid;
}

bool GridForNeighList::GridPosition(XYZ coord, int &ix, int &iy, int &iz) {
    ix = (int)floor(coord[0]/GridSize());
    iy = (int)floor(coord[1]/GridSize());
    iz = (int)floor(coord[2]/GridSize());
    if(ix<0||ix>=Lx()  || iy<0||iy>=Ly()   || iz<0||iz>=Lz())
        return false;
    else
        return true;
}
bool GridForNeighList::AddAtom(AtomInGrid &atom) {
    int ix,iy,iz;
    bool result = GridPosition(atom.xyz,ix,iy,iz);
    if(!result)
        return false;
    else{
        grid[ix][iy][iz].push_back(atom);
        return true;
    }
}
bool GridForNeighList::AddAtom(int index, XYZ xyz, bool ghost,string element) {
    AtomInGrid atom;
    atom.index = index;
    atom.xyz = xyz;
    atom.ghost = ghost;
    atom.element = element;
    return AddAtom(atom);
}
int GridForNeighList::AtomsCount(){
    int count = 0;
    for(int i=0;i<Lx();i++) {
        for (int j = 0; j < Ly(); j++) {
            for (int k = 0; k < Lz(); k++) {
                count += grid[i][j][k].size();
            }
        }
    }
    return count;
}
bool GridForNeighList::NearbyAtoms(XYZ coord, vector<AtomInGrid> &outputVec, bool half) {
    int ix,iy,iz;
    bool result = GridPosition(coord,ix,iy,iz);
    if(!result)
        return false;

    // There are 26 cell plus the central cell, totally 27 cells.
    int cells[27][3];
    bool active[27];
    int counter = 0;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
               cells[counter][0]=i; cells[counter][1]=j; cells[counter][2]=k;
               active[counter]=true;
               counter++;
            }
        }
    }
    // But not all cells are relevant. Firstly, exclude cells on -x/-y/-z directions if half==true
    if(half){
        for(int counter=0;counter<27;counter++){
            int *pCell = cells[counter];
            if(pCell[0]<0 || pCell[1]<0 || pCell[2]<0)
                active[counter] = false;
        }
    }
    // Exclude cells outside the pGrid_
    for(counter=0;counter<27;counter++){
        int *pc = cells[counter];
        pc[0] += ix;
        pc[1] += iy;
        pc[2] += iz;
        if(pc[0]<0 || pc[0]>=Lx() || pc[1]<0 || pc[1]>=Ly() || pc[2]<0 || pc[2] >= Lz())
            active[counter] = false;
    }
    // Finally, report atoms in the active cells
    outputVec.clear();
    for(counter=0;counter<27;counter++){
        if(!active[counter])
            continue;
        int icx = cells[counter][0];
        int icy = cells[counter][1];
        int icz = cells[counter][2];
        auto & vec = grid[icx][icy][icz];
        outputVec.insert(outputVec.end(),vec.begin(),vec.end());
    }
    return true;
}

void GridForNeighList::_showContent_Debug(string dumpfilename) {
    MolecularSystem ms;
    for(int i=0;i<Lx();i++) {
        for (int j = 0; j < Ly(); j++) {
            for (int k = 0; k < Lz(); k++) {
                if(grid[i][j][k].size()>0){
                    cout<<i<<" "<<j<<" "<<k<<":"<<grid[i][j][k].size()<<endl;
                    auto mol = _showGridContentAsMolecule_Debug(i,j,k);
                    ms.AddMolecule(mol);
                }
            }
        }
    }
    ms.RenumberAtomSerials();
    XYZFile xyzfile;
    ms.Write(&xyzfile,dumpfilename);
}

Molecule GridForNeighList::_showGridContentAsMolecule_Debug(int ix, int iy, int iz) {
    vector<AtomInGrid> &vec = grid[ix][iy][iz];
    return _showGridContentAsMolecule_Debug(vec);
}

Molecule GridForNeighList::_showGridContentAsMolecule_Debug(vector<AtomInGrid> &vec) {
    Molecule mol;
    for(int i=0;i<vec.size();i++) {
        Atom a;
        a.element = vec[i].element;
        a.xyz = vec[i].xyz;
        mol.AddAtom(a);
    }
    return mol;
}


double NeighborList::LARGE = 1e9;
double NeighborList::SMALL = 1e-5;

NeighborList::NeighborList(MolecularSystem &ms, double gridSize): pMS_(&ms), gridsize_(gridSize){
    if(gridsize_<0.1)
        error("In NeighborList(), gridsize too small. We suggest 3~5 Angstroms.");
    if(pMS_->Periodic())
        createGridsForPeriodic();
    else
        createGridsForNonPeriodic();
}
void NeighborList::createGridsForNonPeriodic() {
    // Find the minx, miny, minz position of the system.
    // enclose the system with a half-gridsize_ skin,
    // "shift" the system to the new origin_, and create the grids.
    int nAtoms = pMS_->AtomsCount();
    XYZ minXYZ = XYZ(LARGE,LARGE,LARGE);
    XYZ maxXYZ = XYZ(-LARGE,-LARGE,-LARGE);
    for(int i=0;i<nAtoms;i++){
        auto pAtom = pMS_->GetAtom(i);
        for(int j=0;j<3;j++){
            minXYZ[j] = min(minXYZ[j],pAtom->xyz[j]);
            maxXYZ[j] = max(maxXYZ[j],pAtom->xyz[j]);
        }
    }
    // enclose with a skin with thickness = 1/2 of gridsize_
    minXYZ -= XYZ(gridsize_, gridsize_, gridsize_) * 0.5;
    maxXYZ += XYZ(gridsize_, gridsize_, gridsize_) * 0.5;
    systemLengths_ = gridLengths_ = maxXYZ-minXYZ;
    pGrid_ = make_shared<GridForNeighList>(systemLengths_[0], systemLengths_[1], systemLengths_[2], gridsize_);
    gridOrigin_ = minXYZ;
    // Add atoms to the grid
    for(int i=0;i<nAtoms;i++){
        auto pAtom = pMS_->GetAtom(i);
        AtomInGrid a;
        a.index = i;
        a.xyz = pAtom->xyz-minXYZ;
        a.element = pAtom->element;
        a.ghost = false; // no ghost atoms for non-periodic systems
        pGrid_->AddAtom(a);
    }
}
void NeighborList::createGridsForPeriodic() {
    if( ! pMS_->boundary.Orthogonal())
        error("NeighborList::createGridsForPeriodic() only supported for orthogonal system.");
    gridOrigin_ = pMS_->boundary.GetOrigin();
    systemLengths_ = XYZ(pMS_->boundary.GetU()[0],
                             pMS_->boundary.GetV()[1],
                             pMS_->boundary.GetW()[2]);
    if(systemLengths_[0] < SMALL || systemLengths_[1] < SMALL || systemLengths_[2] < SMALL)
        error("In NeighborList::createGridsForPeriodic(), system dimension negative or close to 0.");
    // The size of the grids are wider with skin size = gridsize;
    gridLengths_ = systemLengths_ + XYZ(gridsize_, gridsize_, gridsize_);
    pGrid_ = make_shared<GridForNeighList>(gridLengths_[0], gridLengths_[1], gridLengths_[2], gridsize_);
    // Add atoms to the grid
    int nAtoms = pMS_->AtomsCount();
    for(int i=0;i<nAtoms;i++){
        auto pAtom = pMS_->GetAtom(i);
        AtomInGrid a;
        // atoms must be wrapped in the cell
        a.index = i;
        a.xyz = WrapInCell(pAtom->xyz, GridOrigin(), systemLengths_);
        a.element = pAtom->element;
        a.ghost = false; // first, add the real atom
        pGrid_->AddAtom(a);
        vector<XYZ> equilPos;
        generateEquivalentPositionsForPBC(a.xyz,equilPos);
        for(auto &pos:equilPos){
            pGrid_->AddAtom(i,pos,true,a.element);
        }
    }
}
XYZ NeighborList::WrapInCell(XYZ xyz, XYZ origin, XYZ LxLyLz) {
    xyz = xyz-origin;
    for(int i=0;i<3;i++) {
        while (xyz[i] < 0)
            xyz[i] += LxLyLz[i];
        while (xyz[i] > LxLyLz[i] + SMALL) //+SMALL to void numerical issues: atoms just on the boundary
            xyz[i] -= LxLyLz[i];
    }
    return xyz;
}

void NeighborList::GetNeighbors(XYZ xyz, vector<AtomInGrid> &outputVec, bool half) {
    // Convert original xyz to xyz in the grids. NonPBC and PBC are handled differently.
    if(pMS_->Periodic())
        xyz = WrapInCell(xyz,GridOrigin(),SystemLengths());
    else
        xyz = xyz-GridOrigin();
    pGrid_->NearbyAtoms(xyz,outputVec,half);
    // if the system is NonPBC, or the user needs only a half-neighbor list, that's enough
    if(half || !pMS_->Periodic()) {
        removeAtomListDuplication(outputVec);
        return;
    }

    // The system is PBC and requires full-neighbor list. Need to consider the PBC images
    // if the coord is near the boundary
    vector<XYZ> equilPos;
    generateEquivalentPositionsForPBC(xyz,equilPos);
    for(XYZ &pos:equilPos){
        vector<AtomInGrid> neighbors;
        pGrid_->NearbyAtoms(pos,neighbors,false);
        outputVec.insert(outputVec.end(),neighbors.begin(),neighbors.end());
    }
    // Note the outputVec may contain duplications, now we remove them;
    removeAtomListDuplication(outputVec);
}

double NeighborList::DistanceSquared(XYZ pos1, XYZ pos2) {
    XYZ displ = pos1-pos2;
    if( ! pMS_->Periodic())
        return XYZDot(displ,displ);
    pos1 = WrapInCell(pos1,GridOrigin(),SystemLengths());
    pos2 = WrapInCell(pos2,GridOrigin(),SystemLengths());
    displ = pos1 - pos2;
    for(int i=0;i<3;i++) {
        displ[i] = fabs(displ[i]);
        if(displ[i] > 0.5 * SystemLengths()[i])
            displ[i] = SystemLengths()[i] - displ[i];
    }
    return XYZDot(displ,displ);
}

void NeighborList::generateEquivalentPositionsForPBC(XYZ xyzInCell, vector<XYZ> &outVec) {
    outVec.clear();
    // If the atom is adjacent to the boundary, add its ghost(s)
    if(xyzInCell[0] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], 0, 0));
    if(xyzInCell[1] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(0, systemLengths_[1], 0));
    if(xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(0, 0, systemLengths_[2]));
    // also need the xy, yz, xz, xyz cross terms.
    if(xyzInCell[0] < gridsize_ && xyzInCell[1] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], systemLengths_[1], 0));
    if(xyzInCell[0] < gridsize_ && xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], 0, systemLengths_[2]));
    if(xyzInCell[1] < gridsize_ && xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(0, systemLengths_[1], systemLengths_[2]));
    if(xyzInCell[0] < gridsize_ && xyzInCell[1] < gridsize_ && xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], systemLengths_[1], systemLengths_[2]));
}

void NeighborList::removeAtomListDuplication(vector<AtomInGrid> &vec) {
    vector<AtomInGrid> copy;
    set<int> foundAtoms;
    for(int i=0;i<vec.size();i++) {
        AtomInGrid *pAtom = &vec[i];
        if (foundAtoms.find(vec[i].index) == foundAtoms.end()) {
            foundAtoms.insert(vec[i].index);
            copy.push_back(vec[i]);
        } else {
            continue;
        }
    }
    vec.clear();
    vec.insert(vec.end(),copy.begin(),copy.end());
}

BondDetectorByRules::BondDetectorByRules(double globalCutoff){
    globalCutoff_ = globalCutoff;
    pNlist_ = nullptr;
}
void BondDetectorByRules::ParseFile(string filename){
    ifstream ifs(filename);
    if(!ifs)
        error("Can't open ["+filename+"] to read!");
    string line;
    while(getline(ifs,line)){
        line = StringRemoveComment(line);
        if(line=="")
            continue;
        auto parts = StringSplit(line,',');
        try{
            BondRule rule;
            rule.ele1 = StringToCapitalized(parts[0]);
            rule.ele2 = StringToCapitalized(parts[1]);
            rule.low = stof(parts[2]);
            rule.high = stof(parts[3]);
            rule.type = parts[4];
            rule.low_sqr = pow(rule.low,2);
            rule.high_sqr = pow(rule.high,2);
            rules_.push_back(rule);
        }catch(exception e){
            error("While reading ["+filename+"] at this line:\n"+line);
        }
    }
}
void BondDetectorByRules::ClearRules(){
    rules_.clear();
}
void BondDetectorByRules::Detect(MolecularSystem &ms, bool flushCurrentBonds) {
    // Call this function once to update the date structure
    ms.SearchAtomByGlobalSerial("",true);
    if(flushCurrentBonds)
        FlushAllBonds(ms);

    int nAtoms = ms.AtomsCount();
    pNlist_ = make_shared<NeighborList>(ms,globalCutoff_);

    set<pair<int,int>> foundBonds;
    for(int iFromAtom=0;iFromAtom<nAtoms;iFromAtom++){
        auto pFromAtom = ms.GetAtom(iFromAtom);
        vector<AtomInGrid> candidateAtoms;
        pNlist_->GetNeighbors(pFromAtom->xyz,candidateAtoms,true);
//        if(iFromAtom%100 == 0)
//            cout<<iFromAtom<<endl;
        for(auto &item:candidateAtoms){
            int iToAtom = item.index;
            auto pair1 = make_pair(iFromAtom,iToAtom);
            auto pair2 = make_pair(iToAtom,iFromAtom);
            if(iToAtom == iFromAtom ||
            foundBonds.find(pair1) != foundBonds.end() ||
            foundBonds.find(pair2) != foundBonds.end())
                continue;
            int iBondType = SearchRules(ms,iFromAtom,iToAtom);
            if(iBondType!=-1) {
                AddBond(ms, iFromAtom, iToAtom, this->rules_[iBondType].type);
                foundBonds.insert(pair1);
                foundBonds.insert(pair2);
            }
        }
    }
}
int BondDetectorByRules::SearchRules(MolecularSystem &ms, int iFromAtom, int iToAtom) {
    auto pAtom1 = ms.GetAtom(iFromAtom);
    auto pAtom2 = ms.GetAtom(iToAtom);
    double dist_sqr = pNlist_->DistanceSquared(pAtom1->xyz,pAtom2->xyz);
    int matched_iBond = -1;
    for(int iBond=0;iBond<rules_.size();iBond++){
        BondRule &rule = rules_[iBond];
        if(dist_sqr<rule.low_sqr || dist_sqr>rule.high_sqr)
            continue;
        if( (rule.ele1==pAtom1->element && rule.ele2==pAtom2->element) ||
            (rule.ele1==pAtom2->element && rule.ele2==pAtom1->element) )
            matched_iBond = iBond;
    }
    return matched_iBond;
}
void BondDetectorByRules::FlushAllBonds(MolecularSystem &ms){
    for(int i=0;i<ms.MoleculesCount();i++){
        ms[i].bonds.clear();
    }
    ms.interMolecularBonds.clear();
}
void BondDetectorByRules::AddBond(MolecularSystem &ms,int iFrom, int iTo, std::string bondType) {
    auto pFromAtom = ms.GetAtom(iFrom);
    auto pToAtom = ms.GetAtom(iTo);

    Bond b;
    b.length = sqrt(pNlist_->DistanceSquared(pFromAtom->xyz,pToAtom->xyz));
    b.type = bondType;

    auto iFromMolAtom = ms.SearchAtomByGlobalSerial(pFromAtom->globalSerial);
    auto iToMolAtom = ms.SearchAtomByGlobalSerial(pToAtom->globalSerial);
    if(iFromMolAtom.first == iToMolAtom.first){
        auto mol = ms[iFromMolAtom.first];
        b.atom1 = pFromAtom->serial;
        b.atom2 = pToAtom->serial;
        mol.bonds.push_back(make_shared<Bond>(b));
    }else{
        Bond b;
        b.atom1 = pFromAtom->globalSerial;
        b.atom2 = pToAtom->globalSerial;
        ms.interMolecularBonds.push_back(make_shared<Bond>(b));
    }
}

BondDetectorByRules* GetDefaultBondDetector(double globalCutoff){
    static BondDetectorByRules detector(globalCutoff);
    static bool initialized = false;
    if(!initialized){
        initialized = true;
        detector.ParseFile(DATAFILESPATH+"/bondrules.csv");
    }
    detector.SetGlobalCutoff(globalCutoff);
    return &detector;
}