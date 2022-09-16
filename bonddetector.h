#ifndef CUNIMOLSYS_BONDDETECTOR_H
#define CUNIMOLSYS_BONDDETECTOR_H
#include "universalmolecularsystem.h"

struct AtomInGrid{
    int index; // original index in MolecularSystem
    bool ghost; // whether a ghost atom
    XYZ xyz;   // coordinates in the pGrid_, not in the original system
    string element;
};

class GridForNeighList{
public:
    // Provide the system's dimension in x, y, and z directions
    GridForNeighList(double SystemLx,double SystemLy,double SystemLz,double gridsize);
    ~GridForNeighList();
    // Copy of the grids disallowed.
    GridForNeighList& operator=(const GridForNeighList& o) = delete;
    GridForNeighList(const GridForNeighList& o) = delete;
    inline int Lx(){return Lx_;}
    inline int Ly(){return Ly_;}
    inline int Lz(){return Lz_;}
    inline double GridSize(){return gridsize_;}
    // Calculate the pGrid_ position of coord.
    // if in/out of the pGrid_, returns true/false.
    bool GridPosition(XYZ coord,int &ix,int &iy,int &iz);
    bool AddAtom(AtomInGrid &atom);
    bool AddAtom(int index,XYZ xyz,bool ghost,string element);
    int AtomsCount();
    /* generate the list of atoms near the coord position. this function will return
     * all atoms in the same pGrid_ as well as adjacent grids
     * If half == true, only return grids with larger x,y,and z. */
    bool NearbyAtoms(XYZ coord,vector<AtomInGrid>& output,bool half);
    void _showContent_Debug(string dumpfilename);
    Molecule _showGridContentAsMolecule_Debug(int ix,int iy,int iz);
    Molecule _showGridContentAsMolecule_Debug(vector<AtomInGrid> &vec);
protected: //
    int Lx_;
    int Ly_;
    int Lz_;
    double gridsize_;
    vector<AtomInGrid> ***grid;
};

class NeighborList{
public:
    NeighborList(MolecularSystem &ms,double gridSize);
    inline shared_ptr<GridForNeighList> GetGrid(){return pGrid_;}
    inline XYZ GridOrigin(){return gridOrigin_;}
    inline XYZ GridLengths(){return gridLengths_;}
    inline XYZ SystemLengths(){return systemLengths_;}
    /* Get neighbors. xyz are coordinates in the original system.
     * if half==true, only gives neighbors in the +x, +y, +z direction, useful in neighborlist
     * construction. On the other hand, half==false gives all neighbors, useful in case like
     * clash detection */
    void GetNeighbors(XYZ xyz,vector<AtomInGrid> &outputVec, bool half);
    static XYZ WrapInCell(XYZ xyz,XYZ origin,XYZ LxLyLz);
    // distance^2 of two atoms, possibly in a PBC system.
    double DistanceSquared(XYZ pos1,XYZ pos2);
private:
    /* generate a list of equivalent positions for a coordinate xyzInCell.
     * xyzInCell must be in the cell. If xyzInCell is near the lower boundary of
     * the system, its equivalent pos (excluding the original position)
     * in the x,y,z directions are generated in outVec */
    void generateEquivalentPositionsForPBC(XYZ xyzInCell,vector<XYZ> &outVec);
    void createGridsForPeriodic();
    void createGridsForNonPeriodic();
    void removeAtomListDuplication(vector<AtomInGrid> &vec);
    MolecularSystem *pMS_;
    shared_ptr<GridForNeighList> pGrid_;
    double gridsize_;
    XYZ gridOrigin_;
    XYZ gridLengths_;
    XYZ systemLengths_;
    // For floating point number comparisons
    static double LARGE;
    static double SMALL;
};

struct BondRule{
    string ele1;
    string ele2;
    double low;
    double high;
    string type;
    // automatically calculated
    double low_sqr;
    double high_sqr;
};

class BondDetectorByRules:public BondDetector{
public:
    BondDetectorByRules(double globalCutoff);
    void ClearRules();
    void ParseFile(string filename);
    // Currently flusthCurrentBonds = false not implemented!
    void Detect(MolecularSystem &ms,bool flushCurrentBonds=true) override;
    void Extend(string file_name_of_another_rules_set);
    int SearchRules(MolecularSystem &ms,int iFromAtom, int iToAtom);
    void FlushAllBonds(MolecularSystem &ms);
    void AddBond(MolecularSystem &ms,int iFrom,int iTo, string bondType);
    inline void SetGlobalCutoff(double cutoff) {globalCutoff_ = cutoff;}
private:
    double globalCutoff_;
    vector<BondRule> rules_;
    shared_ptr<NeighborList> pNlist_;
};

BondDetectorByRules* GetDefaultBondDetector(double globalCutoff);

#endif //CUNIMOLSYS_BONDDETECTOR_H
