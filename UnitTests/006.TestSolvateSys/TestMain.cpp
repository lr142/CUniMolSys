#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include "universalmolecularsystem.h"
#include "bonddetector.h"
#include "molecularmanipulator.h"
#include "opener.h"
#include <vector>
#include <fstream>
#include <random>
#include <iomanip>
using namespace std;

TEST(Subsystem,t1){
    MolecularSystem ms;
    QuickOpen(ms,DATAFILESPATH+"/Structures/2-propene.xyz");
    cout<<ms.Summary()<<endl;
    set<string> serials = {"1","2","4"};
    MolSysSubSystemBySerials(ms,serials);
    cout<<ms.Summary()<<endl;
    QuickSave(ms,DATAFILESPATH+"/../dump.mol2");
}

TEST(Solvate,t1){
    MolecularSystem ms;
    QuickOpen(ms,DATAFILESPATH+"/Structures/FuchsSandoff.mol2");

    Boundary region;
    region.SetOrigin({-30,-35,-40});
    XYZ_DTYPE uvw[3][3] = {{50,0,0},{0,50,0},{0,0,40}};
    region.SetUVW(uvw);

    MolSysSolvate(ms,region,WaterType::TIP4P,2.8);
    cout<<ms.Summary()<<endl;

    MolSysReduceToSingleMolecule(ms);
    QuickSave(ms,DATAFILESPATH+"/../dump.mol2");
}
