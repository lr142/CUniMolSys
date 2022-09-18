#ifndef CUNIMOLSYS_MOL2FILE_H
#define CUNIMOLSYS_MOL2FILE_H
#include "universalmolecularsystem.h"
class Mol2File:public MolecularFile{
public:
    Mol2File(bool writeElementInsteadOfType=true);
    bool Read(MolecularSystem &ms, string filename) override;
    bool Write(MolecularSystem &ms, string filename) override;

private:
    bool parseAtomLine(string line,Atom &a);
    bool parseBondLine(string line,Bond &b);
    bool parseBoundaryLine(string line, Boundary &boundary);
    // MOL2文件在输出的时候，只有Name字段与Type字段。很多时候，希望在type字段输出
    // 元素类型，此时可以设置此变量为True。否则正常输出原子类型（例如C.ar）
    bool writeElementInsteadOfType_;
    void showError(int lineno,string line,string filename,string comment,bool fatal);

    void writeAMolecule(Molecule& mol,std::ofstream &ofs);
};
#endif //CUNIMOLSYS_MOL2FILE_H
