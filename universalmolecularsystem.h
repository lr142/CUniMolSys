#ifndef CUNIMOLSYS_UNIVERSALMOLECULARSYSTEM_H
#define CUNIMOLSYS_UNIVERSALMOLECULARSYSTEM_H
#include "utility.h"
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <set>
using std::string;
using std::shared_ptr;
using std::vector;
using std::map;
using std::set;

struct Bond{
    string atom1 = ""; // The (unique) serial of atom1
    string atom2 = ""; // The (unique) serial of atom2
    string type = "";
    double length = 0.0;
};

struct Atom{
    string element = "";
    XYZ xyz = {0,0,0};
    string type = "";
    double charge = 0.0;
    bool flexible = true;
    string name = ""; // User-defined name. May not be unique
    string serial = ""; // ID of the atom within the molecule.
    // serial of each atom must be unique within the molecule.
    // Usually the serial is a number starting from 1 in most
    // molecular file formats.
    string globalSerial = ""; // Must be unique within the
    // whole Molecular System. The MolecularSystem class is
    // responsible for numbering all atoms.
    string layerInfo = ""; // flag for high, low levels in applications like QM/MM
    string parent = ""; // Used in some cases, e.g. the mol it belongs to
    string ShowAsXYZ();
    string ShowAllFields();
};

class Molecule{
public:
    string name = "";
    string serial = ""; // Must be unique within a MolecularSystem
    string type = "";
    string chainID = "";
    string resName = "";
    string resSeq = "";
    // returns ith atom. Supports negative index(counting from end)
    Atom& operator[](int index);
    Bond& GetBond(int index);
    inline int AtomsCount() {return atoms.size();}
    inline int BondsCount() {return bonds.size();}

    inline void AddAtom(const Atom &atom){atoms.push_back(std::make_shared<Atom>(atom));}
    string Summary();

    // Operator = returns a shallow copy of the Molecule
    // Specifically, its atoms and bonds are vectors of pointers pointing to
    // original atoms and bonds
    //Molecule& operator = (const Molecule& other) = default;
    // Returns a "deep" copy of itself.
    Molecule DeepCopy();
public:
    // atoms and bonds, stored as vectors of shared pointers
    // They are implemented as public members for convenience
    // If a Molecule is copied by = operator, its atoms and
    // bonds are not really copied. If "deep" copy is desired,
    // call the Copy() interface.
    vector<shared_ptr<Atom>> atoms;
    vector<shared_ptr<Bond>> bonds;
};

class MolecularSystem;
class Trajectory;

// The abstract base class of molecular file readers
// Concrete molecular file readers are implemented separately.
// The Strategy design pattern is used here.
class MolecularFile{
public:
    virtual bool Read(MolecularSystem &ms, string filename) = 0;
    virtual bool Write(MolecularSystem &ms, string filename) = 0;
    virtual ~MolecularFile(){}
};

// The abstract base class of bond detectors
// Concrete bond detectors are implemented separately.
// The Strategy design pattern is used here.
class BondDetector{
public:
    virtual void Detect(MolecularSystem &ms,bool flushCurrentBonds) = 0;
};

class Boundary{
public:
    Boundary();
    bool Orthogonal();
    void SetUVW(XYZ u,XYZ v,XYZ w);
    void SetUVW(XYZ_DTYPE uvw[3][3]);
    void SetOrigin(XYZ origin);
    inline XYZ GetOrigin(){return origin_;}
    XYZ& operator[](int index);
    void SetLoHi(XYZ_DTYPE lohi[3][2]); // lohi is double[3][2]: {{xlo,xhi},{ylo,yhi},{zlo,zhi}}
    void GetLoHi(XYZ_DTYPE lohi[3][2]);
    inline XYZ GetU(){return uvw[0];}
    inline XYZ GetV(){return uvw[1];}
    inline XYZ GetW(){return uvw[2];}
    void GetUVW(XYZ_DTYPE uvw[3][3]);
    bool Periodic();
    string Show();
private:
    XYZ uvw[3];
    XYZ origin_;
};

class MolecularSystem{
public:
    MolecularSystem(string name="");
    ~MolecularSystem() = default;
    inline void SetName(string n){name=n;}
    inline string GetName(){return name;}
    void Read(MolecularFile *pFile,string filename);
    void Write(MolecularFile *pFile,string filename);
    // Deepcopy the molecularSystem, except the trajactory
    MolecularSystem DeepCopy();
    int AtomsCount();
    int MoleculesCount();
    int BondsCount();
    Molecule& operator[](int index);

    /* Renumber each atom and mol in MolecularSystem.
    * A very important function since continuous and unique numbering are prerequisites to MD and QM calculations.
     * Nearly every molecular reader will call this function after reading the function.
     * Other manipulations on the molecule, such as adding/removing mols/atoms, should also call this function.
    */
    void RenumberAtomSerials(int startingGlobalSerial=1);
    inline void AddMolecule(const Molecule &mol){molecules.push_back(std::make_shared<Molecule>(mol));}
    void DetectBonds(BondDetector *pDetector,bool flushCurrentBonds=true);
    void Clear(); // Discard all molecules/atoms/bonds/trajectories
    void ClearBonds();
    inline bool Periodic(){return boundary.Periodic();}
    string Summary();
    //Geometric Operations:
    void Translate(XYZ offset);
    void Rotate(double clockwise_degree, XYZ axis);
    void FractionalToCartesian();

    /* If this flag is true, the MS has been modified since last time,
     * Adding/removing/reordering atoms/molecules/bonds in the MS will cause this flag to become true (By calling functions
     * such as Clear(), ClearBonds(), AddMolecule(), AddAtom() ).
     * If this flag is true, the user should call RenumberAtomSerials() to reset this flag to false before accessing its
     * contents.  */
    inline bool Dirty(){ return dirty_; }

    /* Below are low level data structures. There are cases (e.g. MolecularFile reader) the caller must work directly
     * with these data structures. However, If the caller modified [molecules] or [interMolecularBonds], it's the
     * caller's responsibility to call RenumberAtomSerials() since in this case the Dirty() flag won't change by itself. */
    Boundary boundary;
    vector<shared_ptr<Molecule>> molecules;
    vector<shared_ptr<Bond>> interMolecularBonds; // 记录分子间的Bond（这些Bond不属于任何一个Molecule）
protected:
    string name;
    bool dirty_; // See comment above;
};

/* This class provides quick access to elements in a MolecularSystem.
 * Because the elements in a MolecularSystem are not indexed, random access, for example search an atom with a
 * given globalSerial, or check whether two atoms in the system are bonded, are slow. In addition, if the system
 * is modified, all indexes become invalid.
 * Therefore, we design this class MolecularSystemIndexer to build indexes of elements in a MolecularSystem
 * and provide quick random access.  but
 * 1. All indexes are built at the creation of this class, the MolecularSystem pass to its constructor should has all
 * atoms and molecules correctly numbers (called MolecularSystem::RenumberAtomSerials)
 * 2. This object becomes invalid when the MolecularSystem has been modified. The user is responsible for not
 * using this object after the MolecularSystem has changed.
 * 3. MolecularSystemAccessor itself will NOT modify the MolecularSystem
 * 4. Most functions in this class would throw exceptions if failed.
 * 4. This class uses the Proxy design pattern (or the Facade pattern?)
 * */
class MolecularSystemAccessor{
public:
    MolecularSystemAccessor(MolecularSystem &ms);
    ~MolecularSystemAccessor() = default;
    Atom& AtomByGlobalIndex(int globalIndex);
    Atom& AtomByMolAndAtomIndex(int molIndex,int atomIndexInMol);
    Atom& AtomByGlobalSerial(string globalSerial);
    int GlobalIndexOfAtom(Atom& atom);
    int GlobalIndexOfAtom(string globalSerial);
    pair<int,int> MolAndLocalIndexOfAtom(Atom& atom);
    pair<int,int> MolAndLocalIndexOfAtom(int globalIndex);
    pair<int,int> MolAndLocalIndexOfAtom(string globalSerial);
    Molecule& MolByIndex(int index);
    Molecule& MolBySerial(string serial);
    int IndexOfMolecule(Molecule& mol);
    Molecule& ParentMolOfAtom(Atom& atom);
    Molecule& ParentMolOfAtom(int globalIndex);
    inline int AtomsCount(){return atoms_count_;}
    inline int MolsCount(){return mols_count_;}
    inline MolecularSystem& GetMolecularSystem(){return ms_;}
    shared_ptr<Bond> GetBond(int iFromAtomIndex,int iToAtomIndex);
    inline bool IfBonded(int iFromAtomIndex,int iToAtomIndex){
        return GetBond(iFromAtomIndex,iToAtomIndex)!=nullptr;
    }
    // SetBond() only affect data within the Accessor, not the actual MolecularSystem!
    void SetBond(int iFromAtomIndex,int iToAtomIndex, Bond &bond);
    inline vector<map<int,shared_ptr<Bond>>>& GetBondedMap(){return bonded_map_;}
private:
    MolecularSystem &ms_;
    int atoms_count_;
    int mols_count_;
    map<string,int> mol_serial_to_index_map_;
    map<string,int> atom_global_serial_to_index_map_;
    vector<int> atom_global_index_to_mol_index_vec_;
    vector<int> atom_global_index_to_local_index_vec_;
    vector<map<string,int>> local_serial_to_local_index_map_;
    /* the function to build the above data structures */
    void indexing_();

    /* bonded_map_ has length of atoms_count_
     * bonded_map_[0] shows the bonding of atom with global index 0; If atom 0 is bonded with atom 2 and 5,
     * then bonded_map_[0] has two items {2, <the pointer to bond 0-2>} and {5, <the pointer to bond 0-5>}
     * Note that this map is symmetric, bond 0-2 will appear both in bonded_map[0] as an item with key=2, and
     * also appear in bonded_map_[2] as an item with key=0. */
    vector<map<int,shared_ptr<Bond>>> bonded_map_;
    // This data structure is built by the following function, which must be called after indexing_()
    void build_bonding_map();
};
#endif //CUNIMOLSYS_UNIVERSALMOLECULARSYSTEM_H
