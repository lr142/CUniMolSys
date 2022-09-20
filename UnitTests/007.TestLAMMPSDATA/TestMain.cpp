#include <iostream>
#include "gtest/gtest.h"
#include "opener.h"
#include "molecularmanipulator.h"
#include "universalmolecularsystem.h"
#include <vector>
#include <fstream>
using namespace std;

class TestLAMMPSData:public ::testing::Test{
protected:
    void SetUp() override{
        string path = DATAFILESPATH+"/../UnitTests/007.TestLAMMPSDATA/";
        QuickOpen(ms[0],path+"relaxed.data");
        //QuickOpen(ms[1],path+"system.data");
    }
    MolecularSystem ms[1];
};

TEST_F(TestLAMMPSData,t1){
    for(int i=0;i<sizeof(ms)/sizeof(ms[0]);i++) {
        int atomCount = ms[i].AtomsCount();
        int molCount = ms[i].MoleculesCount();
        int bondCount = ms[i][0].BondsCount();

        cout<<ms[i].Summary()<<endl;
        QuickSave(ms[i],DATAFILESPATH+"/../dump.mol2");

        ms[i].RenumberAtomSerials(1001);
        EXPECT_EQ(ms[i].AtomsCount(), atomCount);
        EXPECT_EQ(ms[i].MoleculesCount(), molCount);
        EXPECT_EQ(ms[i][0].BondsCount(),bondCount);

    }
//    EXPECT_EQ(ms[0].Periodic(),false);
//    EXPECT_EQ(ms[1].Periodic(),true);
}
