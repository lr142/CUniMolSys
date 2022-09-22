#include "bonddetector.h"
#include "utility.h"
#include "xyzfile.h"
#include "mol2file.h"
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

GridForNeighList::GridForNeighList(double SystemLx, double SystemLy, double SystemLz, double gridsize):gridsize_(gridsize) {
    Lx_ = int(ceil(SystemLx/GridSize()));
    Ly_ = int(ceil(SystemLy/GridSize()));
    Lz_ = int(ceil(SystemLz/GridSize()));
    if(Lx_<=0 or Ly_<=0 or Lz_<=0 or gridsize<=0 ){
        ERROR("In GridForNeighList(), one of the dimension <=0");
    }else if(double(Lx_)*Ly_*Lz_ > 1000000000 ){ // > 1 billion, convert to double to prevent overflow
        ERROR("In GridForNeighList(), too many grids: "+to_string(double(Lx_)*Ly_*Lz_)+", check gridsize_");
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
    if(ix<0 or ix>=Lx() or iy<0 or iy>=Ly() or iz<0 or iz>=Lz())
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
bool GridForNeighList::AddAtom(int index, XYZ xyz, bool ghost) {
    AtomInGrid atom(index,xyz);
    //atom.ghost = ghost;
    //atom.element = element;
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
bool GridForNeighList::NearbyAtoms(XYZ coord, vector<AtomInGrid> &outputVec) {
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
    // But not all cells are relevant.
    // Exclude cells outside the pGrid_
    for(counter=0;counter<27;counter++){
        int *pc = cells[counter];
        pc[0] += ix;
        pc[1] += iy;
        pc[2] += iz;
        if(pc[0]<0 or pc[0]>=Lx() or pc[1]<0 or pc[1]>=Ly() or pc[2]<0 or pc[2] >= Lz())
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

void GridForNeighList::_showContent_Debug(MolecularSystemAccessor &ms, string dumpfilename) {
    MolecularSystem ms_toWrite;
    for(int i=0;i<Lx();i++) {
        for (int j = 0; j < Ly(); j++) {
            for (int k = 0; k < Lz(); k++) {
                if(grid[i][j][k].size()>0){
                    cout<<i<<" "<<j<<" "<<k<<":"<<grid[i][j][k].size()<<endl;
                    auto mol = _debug_ShowGridContentAsMolecule(i, j, k, ms);
                    ms_toWrite.AddMolecule(mol);
                }
            }
        }
    }
    ms_toWrite.RenumberAtomSerials();
    XYZFile xyzfile;
    ms_toWrite.Write(&xyzfile,dumpfilename);
}

Molecule GridForNeighList::_debug_ShowGridContentAsMolecule(int ix, int iy, int iz, MolecularSystemAccessor &ms) {
    vector<AtomInGrid> &vec = grid[ix][iy][iz];
    return _debug_ShowGridContentAsMolecule(vec, ms);
}

Molecule GridForNeighList::_debug_ShowGridContentAsMolecule(vector<AtomInGrid> &vec, MolecularSystemAccessor &ms) {
    Molecule mol;
    for(int i=0;i<vec.size();i++) {
        Atom a;
        a.element = ms.AtomByGlobalIndex(vec[i].index).element;
        a.xyz = vec[i].xyz;
        mol.AddAtom(a);
    }
    return mol;
}

NeighborList::NeighborList(MolecularSystem &ms, double gridSize): msa(ms), gridsize_(gridSize), cached_nlist_(nullptr){
    if(gridsize_<0.1)
        ERROR("In NeighborList(), gridsize too small. We suggest 3~5 Angstroms.");
    if(ms.Periodic())
        createGridsForPeriodic();
    else
        createGridsForNonPeriodic();
}
NeighborList::~NeighborList(){
    ClearCachedNeighborList();
}
void NeighborList::createGridsForNonPeriodic() {
    // Find the minx, miny, minz position of the system.
    // enclose the system with a half-gridsize_ skin,
    // "shift" the system to the new origin_, and create the grids.
    int nAtoms = msa.AtomsCount();
    XYZ minXYZ = XYZ(MY_LARGE,MY_LARGE,MY_LARGE);
    XYZ maxXYZ = XYZ(-MY_LARGE,-MY_LARGE,-MY_LARGE);
    for(int i=0;i<nAtoms;i++){
        Atom &atom = msa.AtomByGlobalIndex(i);
        for(int j=0;j<3;j++){
            minXYZ[j] = min(minXYZ[j],atom.xyz[j]);
            maxXYZ[j] = max(maxXYZ[j],atom.xyz[j]);
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
        Atom &atom = msa.AtomByGlobalIndex(i);
        AtomInGrid a(i,atom.xyz-minXYZ);
        //a.element = pAtom->element;
        //a.ghost = false; // no ghost atoms for non-periodic systems
        pGrid_->AddAtom(a);
    }
}
void NeighborList::createGridsForPeriodic() {
    Boundary & bound = msa.GetMolecularSystem().boundary;
    if( ! bound.Orthogonal())
        ERROR("Only supported for orthogonal system.");
    gridOrigin_ = bound.GetOrigin();
    systemLengths_ = XYZ(bound.GetU()[0],bound.GetV()[1],bound.GetW()[2]);
    if(systemLengths_[0] < MY_SMALL or systemLengths_[1] < MY_SMALL or systemLengths_[2] < MY_SMALL)
        ERROR("System dimension negative or close to 0.");
    // The size of the grids are wider with skin size = gridsize;
    gridLengths_ = systemLengths_ + XYZ(gridsize_, gridsize_, gridsize_);
    pGrid_ = make_shared<GridForNeighList>(gridLengths_[0], gridLengths_[1], gridLengths_[2], gridsize_);
    // Add atoms to the grid
    int nAtoms = msa.AtomsCount();
    for(int i=0;i<nAtoms;i++){
        Atom &atom = msa.AtomByGlobalIndex(i);
        // atoms must be wrapped in the cell
        XYZ pos_in_grid = WrapInCell(atom.xyz, GridOrigin(), systemLengths_);
        AtomInGrid a(i,pos_in_grid);
        //a.element = pAtom->element;
        //a.ghost = false; // first, add the real atom
        pGrid_->AddAtom(a);
        vector<XYZ> equilPos;
        generateEquivalentPositionsForPBC(pos_in_grid,equilPos);
        for(auto &pos:equilPos){
            pGrid_->AddAtom(i,pos,true);
        }
    }
}
XYZ NeighborList::WrapInCell(XYZ xyz, XYZ origin, XYZ LxLyLz) {
    xyz = xyz-origin;
    for(int i=0;i<3;i++) {
        while (xyz[i] < 0)
            xyz[i] += LxLyLz[i];
        while (xyz[i] > LxLyLz[i] + MY_SMALL) //+SMALL to void numerical issues: atoms just on the boundary
            xyz[i] -= LxLyLz[i];
    }
    return xyz;
}

void NeighborList::GetNeighborsFromCoordinates(XYZ xyzInSystem, vector<AtomInGrid> &outputVec) {
    XYZ xyz = SystemToGridXYZ(xyzInSystem);
    pGrid_->NearbyAtoms(xyz,outputVec);
    // if the system is NonPBC, that's enough
    if(!msa.GetMolecularSystem().boundary.Periodic()) {
        removeAtomListDuplication(outputVec);
        return;
    }

    // The system is PBC and requires full-neighbor list. Need to consider the PBC images
    // if the coord is near the boundary
    vector<XYZ> equilPos;
    generateEquivalentPositionsForPBC(xyz,equilPos);
    for(XYZ &pos:equilPos){
        vector<AtomInGrid> neighbors;
        pGrid_->NearbyAtoms(pos,neighbors);
        outputVec.insert(outputVec.end(),neighbors.begin(),neighbors.end());
    }
    // Note the outputVec may contain duplications, now we remove them;
    removeAtomListDuplication(outputVec);
}

void NeighborList::BuildCachedNeighborList(){
    cached_nlist_ = new vector<AtomInGrid>** [pGrid_->Lx()];
    for(int i=0;i<pGrid_->Lx();i++){
        cached_nlist_[i] = new vector<AtomInGrid>* [pGrid_->Ly()];
        for(int j=0;j<pGrid_->Ly();j++){
            cached_nlist_[i][j] = new vector<AtomInGrid> [pGrid_->Lz()];
        }
    }
    for(int i=0;i<pGrid_->Lx();i++){
        for(int j=0;j<pGrid_->Ly();j++){
            for(int k=0;k<pGrid_->Lz();k++){
                XYZ gridPos(i*gridsize_,j*gridsize_,k*gridsize_);
                gridPos += GridOrigin();
                GetNeighborsFromCoordinates(gridPos,cached_nlist_[i][j][k]);
            }
        }
    }
}
void NeighborList::ClearCachedNeighborList() {
    if(cached_nlist_== nullptr)
        return;
    for(int i=0;i<pGrid_->Lx();i++){
        for(int j=0;j<pGrid_->Ly();j++){
            delete [] cached_nlist_[i][j];
        }
        delete [] cached_nlist_[i];
    }
    delete [] cached_nlist_;
    cached_nlist_ = nullptr;
}

vector<AtomInGrid>* NeighborList::GetNeighborsFromCachedLists(XYZ xyzInSystem){
    XYZ xyz = SystemToGridXYZ(xyzInSystem);
    int ix,iy,iz;
    bool result = pGrid_->GridPosition(xyz,ix,iy,iz);
    if(!result || cached_nlist_== nullptr)
        return nullptr;
    else
        return &cached_nlist_[ix][iy][iz];
}

XYZ NeighborList::SystemToGridXYZ(XYZ xyz){
    // Convert original xyz to xyz in the grids. NonPBC and PBC are handled differently.
    if(msa.GetMolecularSystem().boundary.Periodic())
        xyz = WrapInCell(xyz,GridOrigin(),SystemLengths());
    else
        xyz = xyz-GridOrigin();
    return xyz;
}

double NeighborList::DistanceSquared(XYZ pos1, XYZ pos2) {
    XYZ displ = pos1-pos2;
    if( ! msa.GetMolecularSystem().boundary.Periodic())
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
    if(xyzInCell[0] < gridsize_ and xyzInCell[1] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], systemLengths_[1], 0));
    if(xyzInCell[0] < gridsize_ and xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], 0, systemLengths_[2]));
    if(xyzInCell[1] < gridsize_ and xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(0, systemLengths_[1], systemLengths_[2]));
    if(xyzInCell[0] < gridsize_ and xyzInCell[1] < gridsize_ && xyzInCell[2] < gridsize_)
        outVec.push_back(xyzInCell + XYZ(systemLengths_[0], systemLengths_[1], systemLengths_[2]));
}

void NeighborList::removeAtomListDuplication(vector<AtomInGrid> &vec) {
    vector<AtomInGrid> copy;
    set<int> foundAtoms;
    for(int i=0;i<vec.size();i++) {
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

void NeighborList::_debug_ShowNeighbors(XYZ coord,string filename){
    MolecularSystem ms;
    Atom theCoordItself;
    theCoordItself.element = "Pt";
    theCoordItself.xyz = SystemToGridXYZ(coord);
    // Show the neighbors in two ways: directly and from the cached lists:
    // 1. Directly;
    vector<AtomInGrid> vec;
    GetNeighborsFromCoordinates(coord,vec);
    auto mol = pGrid_->_debug_ShowGridContentAsMolecule(vec, msa);
    mol.AddAtom(theCoordItself);
    ms.AddMolecule(mol);

    // 2. From Cached lists:
    auto pVec = GetNeighborsFromCachedLists(coord);
    mol = pGrid_->_debug_ShowGridContentAsMolecule(*pVec, msa);
    mol.AddAtom(theCoordItself);
    ms.AddMolecule(mol);

    XYZFile xyzfile;
    ms.Write(&xyzfile,filename);
}
int NeighborList::_debug_CompareCachedAndDirectNeighborList() {
    int nAtoms = msa.AtomsCount();
    for(int i=0;i<nAtoms;i++) {
        Atom &atom = msa.AtomByGlobalIndex(i);
        vector<AtomInGrid> *fromCached = GetNeighborsFromCachedLists(atom.xyz);
        vector<AtomInGrid> fromDirect;
        GetNeighborsFromCoordinates(atom.xyz, fromDirect);
        if (fromDirect.size() != fromCached->size()) {
            cout << i << "th atom: From Direct: " << fromDirect.size() << "From Cached: " << fromCached->size() << endl;
            return i;
        }
    }
    return -1;
}

BondDetectorByRules::BondDetectorByRules(double globalCutoff){
    SetGlobalCutoff(globalCutoff);
    pNlist_ = nullptr;
}
void BondDetectorByRules::SetGlobalCutoff(double cutoff) {
    globalCutoff_ = cutoff;
    if(cutoff<2.0)
        WARNING("cutoff= "+to_string(cutoff)+" too small, may miss bonds");
    if(cutoff>8.0)
        WARNING("cutoff= " +to_string(cutoff)+" too large, may affect performance.");
}

void BondDetectorByRules::ParseFile(string filename){
    ifstream ifs(filename);
    if(!ifs)
        ERROR("Can't open ["+filename+"] to read!, Please check the DATAFILSPATH in CMakeLists.txt");
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
            string name1 = rule.ele1+"_"+rule.ele2;
            string name2 = rule.ele2+"_"+rule.ele1;
            rules_[name1].push_back(rule);
            if(name1!=name2)
                rules_[name2].push_back(rule);
        }catch(exception e){
            ERROR("While reading ["+filename+"] at this line:\n"+line);
        }
    }
}
void BondDetectorByRules::ClearRules(){
    rules_.clear();
}
void BondDetectorByRules::Detect(MolecularSystem &ms, bool flushCurrentBonds) {
    if(ms.AtomsCount()==0)
        return;
    if(flushCurrentBonds)
        ms.ClearBonds();
    MolecularSystemAccessor msa(ms);
    // msa must be built after ClearBonds()

    int nAtoms = ms.AtomsCount();
    pNlist_ = make_shared<NeighborList>(ms,globalCutoff_);
    // Building Cache
    pNlist_->BuildCachedNeighborList();

    // Debugging, check the consistency, time-consuming:
    /*
    {
        int first_wrong = pNlist_->_debug_CompareCachedAndDirectNeighborList();
        if(first_wrong != -1) {
            // If not equal, check the details for a specific atom
            Atom &atom = msa.AtomByGlobalIndex(first_wrong);
            pNlist_->_debug_ShowNeighbors(atom.xyz, DATAFILESPATH + "/../"
                                                    + atom.globalSerial + ".xyz");
            pNlist_->_debug_CompareCachedAndDirectNeighborList();
            exit(0);
        }
    }*/


    for(int iFromAtom=0;iFromAtom<nAtoms;iFromAtom++){
        Atom &fromAtom = msa.AtomByGlobalIndex(iFromAtom);
        auto candidateAtoms = pNlist_->GetNeighborsFromCachedLists(fromAtom.xyz);

        for(auto &item:*candidateAtoms){
            int iToAtom = item.index;
            Atom &toAtom = msa.AtomByGlobalIndex(iToAtom);
            string elements_names = fromAtom.element + "_" + toAtom.element;
            auto pIter = rules_.find(elements_names);

            if(pIter == rules_.end())
                continue;
            if(iToAtom==iFromAtom or msa.IfBonded(iFromAtom,iToAtom))
                continue;
            int iBondType = searchRules(msa, fromAtom, toAtom, pIter->second);
            if(iBondType!=-1) {
                AddBond(msa, iFromAtom, iToAtom, pIter->second[iBondType].type);
            }
        }
    }
}
int BondDetectorByRules::searchRules(MolecularSystemAccessor &msa, Atom &fromAtom, Atom &toAtom, vector<BondRule> &rule_vec){
    double dist_sqr = pNlist_->DistanceSquared(fromAtom.xyz, toAtom.xyz);
    int matched_iBond = -1;
    for(int iBond=0;iBond<rule_vec.size();iBond++){
        BondRule &rule = rule_vec[iBond];
        if(dist_sqr<rule.low_sqr or dist_sqr>rule.high_sqr)
            continue;
        // Only need to check the distance now
        matched_iBond = iBond;
    }
    return matched_iBond;
}

void BondDetectorByRules::AddBond(MolecularSystemAccessor &msa, int iFromAtom, int iToAtom, string bondType) {
    Atom& fromAtom = msa.AtomByGlobalIndex(iFromAtom);
    Atom& toAtom = msa.AtomByGlobalIndex(iToAtom);
    Bond b;
    b.length = sqrt(pNlist_->DistanceSquared(fromAtom.xyz, toAtom.xyz));
    b.type = bondType;
    int iFromMol = msa.MolAndLocalIndexOfAtom(iFromAtom).first;
    int iToMol = msa.MolAndLocalIndexOfAtom(iToAtom).first;
    if(iFromMol == iToMol){
        b.atom1 = fromAtom.serial;
        b.atom2 = toAtom.serial;
        msa.GetMolecularSystem()[iFromMol].bonds.push_back(make_shared<Bond>(b));
    }else{
        Bond b;
        b.atom1 = fromAtom.globalSerial;
        b.atom2 = toAtom.globalSerial;
        msa.GetMolecularSystem().interMolecularBonds.push_back(make_shared<Bond>(b));
    }
    msa.SetBond(iFromAtom,iToAtom,b);
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
