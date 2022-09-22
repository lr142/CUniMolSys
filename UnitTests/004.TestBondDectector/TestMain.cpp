#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include "mol2file.h"
#include "xyzfile.h"
#include "universalmolecularsystem.h"
#include "bonddetector.h"
#include <vector>
#include <fstream>
using namespace std;

class TestSplitting: public ::testing::Test{
protected:
    void SetUp() override{
        string paths[] = {
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
    }
    vector<MolecularSystem> ms;
    Mol2File mol2file;
    XYZFile xyzfile;

    MolecularSystem SelectRegionFromMolecularSystem(MolecularSystem &ms, XYZ origin, XYZ lengths){
        MolecularSystem newMS;
        newMS.boundary.SetOrigin(origin);
        XYZ_DTYPE uvw[3][3] = {{lengths[0],0,0},{0,lengths[1],0},{0,0,lengths[2]}};
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

// Test if there are memory leak
//TEST_F(TestSplitting,memory){
//    while(true) {
//        GridForNeighList pGrid_(50, 50, 50, 0.37);
//        int nAtom = ms[1].AtomsCount();
//        for (int i = 0; i < nAtom; i++) {
//            auto pAtom = ms[1].GetAtom(i);
//            pGrid_.AddAtom(i, pAtom->xyz, false, pAtom->element);
//        }
//    }
//}

// Test if the GridForNeighList works fine
//TEST_F(TestSplitting,t1){
//    GridForNeighList pGrid_(50,50,50,1.4);
//    int nAtom = ms[1].AtomsCount();
//    for(int i=0;i<nAtom;i++){
//        auto pAtom = ms[1].GetAtom(i);
//        bool result = pGrid_.AddAtom(i,pAtom->xyz,false,pAtom->element);
//    }
//    pGrid_._showContent_Debug(DATAFILESPATH+"/../dump.xyz");
//}

// Test if the GridForNeighList works fine
//TEST_F(TestSplitting,t2){
//    GridForNeighList pGrid_(50,50,50,4.5);
//    int nAtom = ms[1].AtomsCount();
//    for(int i=0;i<nAtom;i++){
//        auto pAtom = ms[1].GetAtom(i);
//        bool result = pGrid_.AddAtom(i,pAtom->xyz,false,pAtom->element);
//    }
//    MolecularSystem newMS;
//    double meshgrid[3] = {4.3,4.71,4.46};
//    for(double x=0; x < 50; x+=meshgrid[0]){
//        cout<<"x="<<x<<endl;
//        for(double y=0; y < 50; y+=meshgrid[1]){
//            //cout<<"y="<<y<<endl;
//            for(double z=0; z < 50; z+=meshgrid[2]){
//                vector<AtomInGrid> vec;
//                pGrid_.NearbyAtoms(XYZ(x, y, z), vec, false);
//                Molecule mol = pGrid_._debug_ShowGridContentAsMolecule(vec);
//                newMS.AddMolecule(mol);
//            }
//        }
//    }
//    newMS.RenumberAtomSerials();
//    cout<<newMS.Summary()<<endl;
//    newMS.Write(&xyzfile,DATAFILESPATH+"/../dump.xyz");
//}

// Test the Neighborlist generation
//TEST_F(TestSplitting,neighlist1){
//    MolecularSystem thems =
//            SelectRegionFromMolecularSystem(ms[1],
//            XYZ(0,0,0),XYZ(50,50,50));
//    cout<<thems.Summary()<<endl;
//
//    NeighborList nlist(thems,4.5);
//    int nAtom = thems.AtomsCount();
//    MolecularSystem toWrite;
//    for(double x=0;x<50;x+=6){
//        for(double y=0;y<50;y+=6){
//            for(double z=0;z<50;z+=6){
//                vector<AtomInGrid> outlist;
//                nlist.GetNeighborsFromCoordinates(XYZ(x,y,z),outlist,true);
//                auto mol = nlist.GetGrid()->_debug_ShowGridContentAsMolecule(outlist);
//                toWrite.AddMolecule(mol);
//            }
//        }
//    }
//    toWrite.Write(&xyzfile,DATAFILESPATH+"/../dump.xyz");
//    cout<<toWrite.Summary()<<endl;
//}

// Test Bond Detection
TEST_F(TestSplitting, findbonds){
    for(int iNo=0;iNo<ms.size();iNo++) {
        MolecularSystem &thems = ms[iNo];
        thems.RenumberAtomSerials();
        cout<<thems.Summary()<<endl;
        int expected_bond_count = thems.BondsCount();
        for (double cutoff = 2; cutoff < 8; cutoff += 0.5) {
            auto start = chrono::system_clock::now();
            thems.DetectBonds(GetDefaultBondDetector(cutoff), true);
            EXPECT_EQ(thems.BondsCount(), expected_bond_count);
            auto end = chrono::system_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            cout << "Cutoff = " << cutoff << " consumes " << duration.count() << "ms" << endl;
        }

        int bondsCountToDelete = thems.BondsCount()/2;
        thems.molecules[0]->bonds.erase(thems.molecules[0]->bonds.begin() ,
                                        thems.molecules[0]->bonds.begin() + bondsCountToDelete);
        EXPECT_EQ(thems.BondsCount(), expected_bond_count - bondsCountToDelete);
        thems.DetectBonds(GetDefaultBondDetector(4.), false);
        EXPECT_EQ(thems.BondsCount(), expected_bond_count);

        thems.Write(&mol2file, DATAFILESPATH + "/../dump" + to_string(iNo)+".mol2");
    }
}