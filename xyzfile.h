#ifndef CUNIMOLSYS_XYZFILE_H
#define CUNIMOLSYS_XYZFILE_H
#include "universalmolecularsystem.h"
class XYZFile:public MolecularFile{
public:
    bool Read(MolecularSystem &ms, string filename) override;
    bool Write(MolecularSystem &ms, string filename) override;
};
#endif //CUNIMOLSYS_XYZFILE_H
