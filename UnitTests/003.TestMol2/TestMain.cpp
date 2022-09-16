#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include "mol2file.h"
#include "universalmolecularsystem.h"
#include <vector>
#include <fstream>
using namespace std;

class TestMol2:public ::testing::Test{
protected:
    void SetUp() override{
        ms[0].Read(&mol2file,DATAFILESPATH+"/Structures/FuchsSandoff.mol2");
        cout<<ms[0].Summary()<<endl;
        ms[1].Read(&mol2file,DATAFILESPATH+"/Structures/tip4p_water_5nm_chunk.mol2");
        cout<<ms[1].Summary()<<endl;
    }
    MolecularSystem ms[2];
    Mol2File mol2file;
};

TEST_F(TestMol2,t1){
    for(int i=0;i<2;i++) {
        int atomCount = ms[i].AtomsCount();
        int molCount = ms[i].MoleculesCount();
        int bondCount = ms[i][0].BondsCount();

        ms[i].RenumberAtomSerials(1001);
        EXPECT_EQ(ms[i].AtomsCount(), atomCount);
        EXPECT_EQ(ms[i].MoleculesCount(), molCount);
        EXPECT_EQ(ms[i][0].BondsCount(),bondCount);


        string toWrite = DATAFILESPATH+"/../dump"+to_string(i+1)+".mol2";
        ms[i].Write(&mol2file,toWrite);
        ms[i].Read(&mol2file,toWrite);

        EXPECT_EQ(ms[i].AtomsCount(), atomCount);
        EXPECT_EQ(ms[i].MoleculesCount(), molCount);
        EXPECT_EQ(ms[i][0].BondsCount(),bondCount);

        ms[i].RenumberAtomSerials(1);
        ms[i].Write(&mol2file,toWrite);

    }
    EXPECT_EQ(ms[0].Periodic(),false);
    EXPECT_EQ(ms[1].Periodic(),true);
}
