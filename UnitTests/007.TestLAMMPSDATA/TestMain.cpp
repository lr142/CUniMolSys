#include <iostream>
#include "gtest/gtest.h"
#include "opener.h"
#include "universalmolecularsystem.h"
#include <vector>
#include <fstream>
using namespace std;

class TestLAMMPSData:public ::testing::Test{
protected:
    void SetUp() override{
        string path = DATAFILESPATH+"/../UnitTests/007.TestLAMMPSDATA/";
        QuickOpen(ms[0],path+"test.data");
        cout<<ms[0].Summary()<<endl;
        QuickOpen(ms[1],path+"system.data");
    }
    MolecularSystem ms[2];
};

TEST_F(TestLAMMPSData,t1){
    for(int i=0;i<2;i++) {
        int atomCount = ms[i].AtomsCount();
        int molCount = ms[i].MoleculesCount();
        int bondCount = ms[i][0].BondsCount();

        ms[i].RenumberAtomSerials(1001);
        EXPECT_EQ(ms[i].AtomsCount(), atomCount);
        EXPECT_EQ(ms[i].MoleculesCount(), molCount);
        EXPECT_EQ(ms[i][0].BondsCount(),bondCount);

    }
//    EXPECT_EQ(ms[0].Periodic(),false);
//    EXPECT_EQ(ms[1].Periodic(),true);
}
