#ifndef CUNIMOLSYS_BONDDETECTOR_H
#define CUNIMOLSYS_BONDDETECTOR_H
#include "universalmolecularsystem.h"

struct AtomInGrid{
    int index; // original index in MolecularSystem
    //bool ghost; // whether a ghost atom
    XYZ xyz;   // coordinates in the pGrid_, not in the original system
    //string element;
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
    bool AddAtom(int index,XYZ xyz,bool ghost);
    int AtomsCount();
    /* generate the list of atoms near the coord position. this function will return
     * all atoms in the same pGrid_ as well as adjacent grids */
    bool NearbyAtoms(XYZ coord,vector<AtomInGrid>& output);
    void _showContent_Debug(MolecularSystemAccessor &ms, string dumpfilename);
    Molecule _debug_ShowGridContentAsMolecule(int ix, int iy, int iz, MolecularSystemAccessor &ms);
    Molecule _debug_ShowGridContentAsMolecule(vector<AtomInGrid> &vec, MolecularSystemAccessor &ms);
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
    ~NeighborList();
    NeighborList(const NeighborList& nlist) = delete;
    NeighborList& operator=(const NeighborList& nlist) = delete;
    inline shared_ptr<GridForNeighList> GetGrid(){return pGrid_;}
    inline XYZ GridOrigin(){return gridOrigin_;}
    inline XYZ GridLengths(){return gridLengths_;}
    inline XYZ SystemLengths(){return systemLengths_;}
    /* Get neighbors. xyz are coordinates in the original system.
     * Note: this per-atom operation is time-consuming */
    void GetNeighborsFromCoordinates(XYZ xyzInSystem, vector<AtomInGrid> &outputVec);

    /* Assume all atoms in the same grid have the same neighbor list
     * To speed-up the process of fetching neighbor list, the BuildCache function is called. */
    void BuildCachedNeighborList();
    void ClearCachedNeighborList();
    vector<AtomInGrid>* GetNeighborsFromCachedLists(XYZ xyzInSystem);

    XYZ SystemToGridXYZ(XYZ xyzInSys);
    static XYZ WrapInCell(XYZ xyzInSystem,XYZ origin,XYZ LxLyLz);
    // distance^2 of two atoms, possibly in a PBC system.
    double DistanceSquared(XYZ pos1,XYZ pos2);
    void _debug_ShowNeighbors(XYZ coord, string filename);
    int _debug_CompareCachedAndDirectNeighborList();
private:
    /* generate a list of equivalent positions for a coordinate xyzInCell.
     * xyzInCell must be in the cell. If xyzInCell is near the lower boundary of
     * the system, its equivalent pos (excluding the original position)
     * in the x,y,z directions are generated in outVec */
    void generateEquivalentPositionsForPBC(XYZ xyzInCell,vector<XYZ> &outVec);
    void createGridsForPeriodic();
    void createGridsForNonPeriodic();
    void removeAtomListDuplication(vector<AtomInGrid> &vec);
    MolecularSystemAccessor msa;
    shared_ptr<GridForNeighList> pGrid_;
    double gridsize_;
    XYZ gridOrigin_;
    XYZ gridLengths_;
    XYZ systemLengths_;
    vector<AtomInGrid> *** cached_nlist_;
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
    void Detect(MolecularSystem &ms,bool flushCurrentBonds=true) override;
    void AddBond(MolecularSystemAccessor &msa, int iFromAtom, int iToAtom, string bondType);
    void SetGlobalCutoff(double cutoff);
private:
    int searchRules(MolecularSystemAccessor &msa, Atom &fromAtom, Atom &toAtom, vector<BondRule> &rule_vec);
    double globalCutoff_;
    map<string,vector<BondRule>> rules_;
    shared_ptr<NeighborList> pNlist_;
};

BondDetectorByRules* GetDefaultBondDetector(double globalCutoff);

#endif //CUNIMOLSYS_BONDDETECTOR_H
