#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include "mol2file.h"
#include "vmdinterface.h"
#include "opener.h"
#include "universalmolecularsystem.h"
#include <vector>
#include <fstream>
using namespace std;

class TestVMD:public ::testing::Test{
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

TEST_F(TestVMD,t1){
    bool periodic[2] = {false,true}; // sys0 non-pbc, sys1 pbc
    for(int i=0;i<2;i++) {

        string no = to_string(i);
        QuickSave(ms[i], no + ".mol2");
        VMDInterface vi(no + ".mol2", no + "vmd.tcl", false);

        EXPECT_EQ(ms[i].Periodic(), periodic[i]);
        if(periodic[i]) {
            //ms[i].boundary.SetOrigin({-2,-3.5,-4});
            vi.AddPeriodicBoundaryBox(ms[i].boundary, "blue", 2, "dashed");
        }else {
            int nAtoms = ms[i].AtomsCount();
            vector<int> indexes(nAtoms);
            vector<string> labels(nAtoms);
            MolecularSystemAccessor msa(ms[i]);
            for (int iAtom = 0; iAtom < nAtoms; iAtom++) {
                indexes[iAtom] = iAtom;
                labels[iAtom] = msa.AtomByGlobalIndex(iAtom).element;
            }
            double offset[3] = {0.06, 0.07, 0.08};
            vi.AddLabels(indexes, labels, offset);
        }

        vi.Run();
    }

}
