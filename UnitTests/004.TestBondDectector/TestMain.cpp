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

class TestBD:public ::testing::Test{
protected:
    void SetUp() override{
        ms[0].Read(&mol2file,DATAFILESPATH+"/Structures/FuchsSandoff.mol2");
        cout<<ms[0].Summary()<<endl;
        ms[1].Read(&mol2file,DATAFILESPATH+"/Structures/tip4p_water_5nm_chunk.mol2");
        cout<<ms[1].Summary()<<endl;
    }
    MolecularSystem ms[2];
    Mol2File mol2file;
    XYZFile xyzfile;

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

// Test if there are memory leak
//TEST_F(TestBD,memory){
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
//TEST_F(TestBD,t1){
//    GridForNeighList pGrid_(50,50,50,1.4);
//    int nAtom = ms[1].AtomsCount();
//    for(int i=0;i<nAtom;i++){
//        auto pAtom = ms[1].GetAtom(i);
//        bool result = pGrid_.AddAtom(i,pAtom->xyz,false,pAtom->element);
//    }
//    pGrid_._showContent_Debug(DATAFILESPATH+"/../dump.xyz");
//}

// Test if the GridForNeighList works fine
//TEST_F(TestBD,t2){
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
//                Molecule mol = pGrid_._showGridContentAsMolecule_Debug(vec);
//                newMS.AddMolecule(mol);
//            }
//        }
//    }
//    newMS.RenumberAtomSerials();
//    cout<<newMS.Summary()<<endl;
//    newMS.Write(&xyzfile,DATAFILESPATH+"/../dump.xyz");
//}

// Test the Neighborlist generation
//TEST_F(TestBD,neighlist1){
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
//                nlist.GetNeighbors(XYZ(x,y,z),outlist,true);
//                auto mol = nlist.GetGrid()->_showGridContentAsMolecule_Debug(outlist);
//                toWrite.AddMolecule(mol);
//            }
//        }
//    }
//    toWrite.Write(&xyzfile,DATAFILESPATH+"/../dump.xyz");
//    cout<<toWrite.Summary()<<endl;
//}

// Test Bond Detection
TEST_F(TestBD,findbonds){
    MolecularSystem thems =
            SelectRegionFromMolecularSystem(ms[1],
            XYZ(0,0,0),XYZ(50,50,50));
    thems.RenumberAtomSerials();
    thems.Summary();
    thems.DetectBonds(GetDefaultBondDetector(5.0));
    thems.Summary();
}