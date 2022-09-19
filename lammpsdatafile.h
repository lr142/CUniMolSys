#ifndef CUNIMOLSYS_LAMMPSDATAFILE_H
#define CUNIMOLSYS_LAMMPSDATAFILE_H
#include "universalmolecularsystem.h"

class LAMMPSDataFile:public MolecularFile{
public:
    LAMMPSDataFile() = default;
    bool Read(MolecularSystem &ms, string filename) override;
    bool Write(MolecularSystem &ms, string filename) override;
};

#endif //CUNIMOLSYS_LAMMPSDATAFILE_H
