#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include "mol2file.h"
#include "xyzfile.h"
#include "universalmolecularsystem.h"
#include "bonddetector.h"
#include "molecularmanipulator.h"
#include <vector>
#include <fstream>
#include <random>
using namespace std;

class TestBD:public ::testing::Test{
protected:
    void SetUp() override{
        string paths[] = {
                DATAFILESPATH+"/Structures/Given1.mol2", //Special Empty Case. See Below.
                DATAFILESPATH+"/Structures/Given1.mol2",
                DATAFILESPATH+"/Structures/Given2.mol2",
                DATAFILESPATH+"/Structures/FuchsSandoff.mol2",
                DATAFILESPATH+"/Structures/PolyPS.mol2",
                DATAFILESPATH+"/Structures/tip4p_water_5nm_chunk.mol2",
                DATAFILESPATH+"/Structures/spc_water_5nm_chunk.mol2",
        };
        for(auto path:paths) {
            ms.push_back(MolecularSystem());
            ms[ms.size() - 1].Read(&mol2file, path);
        }
        // ms[0] is a special case: empty system
        ms[0].Clear();ms[0].ClearBonds();

        e.seed(0);
    }
    vector<MolecularSystem> ms;
    Mol2File mol2file;
    XYZFile xyzfile;
    default_random_engine e;

    MolecularSystem SelectRegionFromMolecularSystem(MolecularSystem &ms, XYZ origin, XYZ lengths){
        MolecularSystem newMS;
        newMS.boundary.SetOrigin(origin);
        double uvw[3][3] = {{lengths[0],0,0},{0,lengths[1],0},{0,0,lengths[2]}};
        newMS.boundary.SetUVW(uvw);

        for(int i=0;i<ms.MoleculesCount();i++){
            Molecule newMol;
            auto curMol = ms[i];
            for(int j=0;j<curMol.AtomsCount();j++){
                auto curAtom = curMol[j];
                XYZ pos = curAtom.xyz - origin;
                if(pos[0]<lengths[0] && pos[1]<lengths[1] && pos[2]<lengths[2])
                    newMol.AddAtom(curAtom);
            }
            newMS.AddMolecule(newMol);
        }
        newMS.RenumberAtomSerials();
        return newMS;
    }

};

// Test Bond Detection
TEST_F(TestBD,findbonds){

    uniform_real_distribution<double> urd(-5,5);

    for(int iNo=0;iNo<ms.size();iNo++) {

        cout<<"Test Case "<<iNo+1<<" of "<<ms.size()<<":"<<endl;

        MolecularSystem &thems = ms[iNo];
        thems.RenumberAtomSerials();
        cout<<"Original: "<<thems.Summary()<<endl;
        int nAtoms = thems.AtomsCount();
        int nBonds = thems.BondsCount();

        vector<int> pieces;
        for(int n=1;n<thems.AtomsCount()*1.5;n+=thems.AtomsCount()/7){
            pieces.push_back(n);
        }
        pieces.push_back(1);

        for(int n:pieces){
            MolSysRandomSplit(thems, n);
            // Randomly split
            EXPECT_EQ(thems.AtomsCount(),nAtoms);
            EXPECT_EQ(thems.BondsCount(),nBonds);
            EXPECT_EQ(thems.MoleculesCount(),nAtoms==0?0:n);
            thems.DetectBonds(GetDefaultBondDetector(5.0),true);
            cout << "Splitted: \t"<<thems.Summary() << endl;
            // Then auto detect
            EXPECT_EQ(thems.BondsCount(),nBonds);
            cout << "AutoDetect: \t"<<thems.Summary() << endl;

            // Random Translate
            thems.Translate(XYZ(urd(e),urd(e),urd(e)));
            thems.DetectBonds(GetDefaultBondDetector(5.0),true);
            EXPECT_EQ(thems.BondsCount(),nBonds);
            cout << "Shifted: \t"<<thems.Summary() << endl;

        }

        thems.Write(&mol2file, DATAFILESPATH + "/../dump" + to_string(iNo)+".mol2");
    }
}