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
    // 一个辅助函数，返回一个连结列表(a list of lists)，例如
    // GetBondedMap[0] 记录与index（非serial）=0 的原子连接的原子的indexes
    // 如果有悬空键（即有的键只有一端，没有另一端）,则出错。
    // 此数据结构较为耗费资源，故只有传递 update = true时，才生成，
    // 若update = false,仅返回上次计算的结果（以static形式存储）
    vector<vector<int>>& GetBondedMap(bool update);
    inline void AddAtom(const Atom &atom){atoms.push_back(std::make_shared<Atom>(atom));}
    string Summary();

    // Operator = returns a shallow copy of the Molecule
    // Specifically, its atoms and bonds are vectors of pointers pointing to
    // original atoms and bonds
    //Molecule& operator = (const Molecule& other) = default;
    // Returns a "deep" copy of itself.
    Molecule DeepCopy();

    /* 分子内部结构合理性（一致性）检测。如果所有检查都通过，则返回 True，否则返回 False
    # 检查的内容包含以下几项：
    # 1. Atom序列号必须唯一
    # 2. 没有悬空键。通过使用辅助函数BondedMap来实现*/
    bool ConsistencyCheck();
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
};

// The abstract base class of bond detectors
// Concrete bond detectors are implemented separately.
// The Strategy design pattern is used here.
class BondDetector{
public:
    virtual void Detect(MolecularSystem &ms,bool flushCurrentBonds) = 0;
};

class Boundary{
/* # 边界条件。非周期性体系可设为None。对周期性系，可以使用以下两种惯例之一：
# 1. 设为3x3的矩阵，代表u,v,w三个lattice向量。此惯例可以处理非Orthogonal的晶胞, 记录于boundaryInUVW
# 2. 对orthogonal体系，特别是用于LAMMPS时，使用3x2矩阵，代表[[xlo,xhigh],[ylo,yhigh],[zlo,zhigh]]. 这种方式也是LAMMPSDATAFile
# 读取文件时采用的表示方式。记录于 boundaryInLoHi。原点记录于 origin_
# 具体使用哪种惯例依调用者的需要而定。
# 对orthogonal 体系， origin_ + boundaryInUVW <--> boundaryInLoHi可以互转 */
public:
    Boundary(bool orthogonal = true);
    inline bool Orthogonal(){return orthogonal_;}
    void SetUVW(XYZ uvw[3]);
    void SetUVW(XYZ u,XYZ v,XYZ w);
    void SetUVW(double uvw[3][3]);
    void SetOrigin(XYZ origin);
    inline XYZ GetOrigin(){return origin_;}
    XYZ& operator[](int index);
    void SetLoHi(double lohi[3][2]); // lohi is double[3][2]: {{xlo,xhi},{ylo,yhi},{zlo,zhi}}
    void GetLoHi(double lohi[3][2]);
    inline XYZ GetU(){return uvw[0];}
    inline XYZ GetV(){return uvw[1];}
    inline XYZ GetW(){return uvw[2];}
    void GetUVW(double uvw[3][3]);
    bool Periodic();
    string Show();
private:
    bool orthogonal_;
    XYZ uvw[3];
    XYZ origin_;
};

class MolecularSystem{
public:
    MolecularSystem(string name="");
    ~MolecularSystem() = default;
    void Read(MolecularFile *pFile,string filename);
    void Write(MolecularFile *pFile,string filename);
    // Deepcopy the molecularsystem, including the trajectory
    MolecularSystem DeepCopyWithTrajectory();
    // Deepcopy the molecularSystem, except the trajactory
    MolecularSystem DeepCopy();
    int AtomsCount();
    int MoleculesCount();
    int BondsCount();
    Molecule& operator[](int index);
    /* Get Atom by index in the whole molecular system. Supports negative Index.
    * Returns nullptr if the index is out of range*/
    shared_ptr<Atom> GetAtom(int index);
    /* According to molecule's serial, return molecule's index.
     * to speedup, the function builds cached data structures upon first call
     * the cached data structure will be invalid if the MolSys has been modified
     * The user should set update=true in those situations to update the cached
     * data. Returns <-1,-1> if no such serial exists
     * # Warning: Call with update=false only when you're absolutely sure!  */
    int SearchMoleculeBySerial(string molSerial,bool update=true);
    /* According to an atom's globalSerial, return the index of its molecule.
     * and the index of the atom within this molecule.
     * to speedup, the function builds cached data structures upon first call
     * the cached data structure will be invalid if the MolSys has been modified
     * The user should set update=true in those situations to update the cached
     * data. Returns <-1,-1> if no such serial exists
     * # Warning: Call with update=false only when you're absolutely sure!
     * */
    std::pair<int,int> SearchAtomByGlobalSerial(string atomGlobalSerial,bool update=true);
    /*  # 重新编号 MolecularSystem中每个原子的serials 和 globalSerial。
    # 这是一个非常重要的函数，因为很多情况下，原子编号的唯一性与连续性是开展MD和QM计算的必要条件。
    # 原则上，几乎所有的MolecularFile(文件读写器)在读取一个MolecularSystem结束时，都会自发调用此函数来对原子统一编号，只有少数例外，
    # 例如LAMMPSDATAFile，因为此时文件中已经有原子统一编号,可能需要保留，不再重新编号。
    # 在对MolecularSystem作其它操作，例如分解、合并、添加成份等后，一般也需要调用此函数。
    # Warning: Call with update=false only when you're absolutely sure!
    */
    inline void AddMolecule(const Molecule &mol){molecules.push_back(std::make_shared<Molecule>(mol));}
    void RenumberAtomSerials(int startingGlobalSerial=1);
    void DetectBonds(BondDetector *pDetector,bool flushCurrentBonds=true);
    void Clear(); // Discard all molecules/atoms/bonds/trajectories
    inline bool Periodic(){return boundary.Periodic();}
    string Summary();
    //Geometric Operations:
    void Translate(XYZ offset);
    void Rotate(double clockwise_degree, XYZ axis);
    void FractionalToCartesian();

public:
    string name = "";
    Boundary boundary;
    vector<shared_ptr<Molecule>> molecules;
    vector<shared_ptr<Bond>> interMolecularBonds; // 记录分子间的Bond（这些Bond不属于任何一个Molecule）
    // # 体系的演化轨迹（如有）
    shared_ptr<Trajectory> trajectory = nullptr;
};
#endif //CUNIMOLSYS_UNIVERSALMOLECULARSYSTEM_H
