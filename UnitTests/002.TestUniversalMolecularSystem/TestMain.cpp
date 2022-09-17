#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include "xyzfile.h"
#include "universalmolecularsystem.h"
#include <vector>
#include <fstream>
using namespace std;

//TEST(TestBoundary,t1){
//    Boundary b(true);
//    double lohi[3][2] = {{-5,10},{0,20},{100,200}};
//    b.SetLoHi(lohi);
//    cout<<b.Show()<<endl;
//    cout<<b[0]<<endl;
//    b.SetUVW(b[0]+XYZ(1,0,0),b[1]+XYZ(0,2,0),b[2]+XYZ(0,0,3));
//    cout<<b.Show()<<endl;
//    b.SetOrigin(XYZ(0,0,0));
//    cout<<b.Show()<<endl;
//}

class TestUniversalMolecularSystem:public ::testing::Test{
protected:
    void SetUp() override{
        ms.Read(&xyzfile,DATAFILESPATH+"/Structures/FuchsSandoff.xyz");
        cout<<ms.Summary()<<endl;
        ms2.Read(&xyzfile,DATAFILESPATH+"/Structures/2-propene.xyz");
    }
    MolecularSystem ms;
    MolecularSystem ms2;
    XYZFile xyzfile;
};

TEST_F(TestUniversalMolecularSystem,t1){
    MolecularSystem ms_composite = ms.DeepCopy();
    XYZ axis(1,1,2);
    for (int i = 0; i < 360; i += 10) {
        MolecularSystem ms_copy = ms.DeepCopy();
        ms_copy.Rotate(i, axis);
        ms_composite.molecules.push_back(ms_copy.molecules[0]);
        ms_composite.RenumberAtomSerials();
    }
    string dumpfile = DATAFILESPATH+"/../dump.xyz";
    ms_composite.Write(&xyzfile,dumpfile);
    ms_composite.Read(&xyzfile,dumpfile);
    cout<<ms_composite.Summary()<<endl;
    EXPECT_EQ(ms_composite.AtomsCount(),9250);
}

TEST_F(TestUniversalMolecularSystem,t2){
    MolecularSystem molSys = ms.DeepCopy();
    molSys.RenumberAtomSerials();
    molSys.AddMolecule(ms2.molecules[0]->DeepCopy());
    molSys.RenumberAtomSerials();
    MolecularSystemAccessor msa(molSys);
    cout<<molSys.Summary()<<endl;
    for(int i=0;i<molSys.AtomsCount();i++){
        auto index = msa.MolAndLocalIndexOfAtom(to_string(i+1));
        EXPECT_EQ(
                molSys[index.first][index.second].globalSerial,
                msa.AtomByGlobalIndex(i).globalSerial
        );
    }
}