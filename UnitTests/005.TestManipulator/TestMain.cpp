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

class TestSplitting: public ::testing::Test{
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

        e.seed(1);
    }
    vector<MolecularSystem> ms;
    Mol2File mol2file;
    XYZFile xyzfile;
    default_random_engine e;

    MolecularSystem SelectRegionFromMolecularSystem(MolecularSystem &ms, XYZ origin, XYZ lengths);
    void ConnectEveryAtomInEachMolecule(MolecularSystem &ms);
};
// Select a specific region in the MolSys.
MolecularSystem TestSplitting::SelectRegionFromMolecularSystem(MolecularSystem &ms, XYZ origin, XYZ lengths){
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
// Artificially connect all atoms within a molecule
void TestSplitting::ConnectEveryAtomInEachMolecule(MolecularSystem &ms){
    MolecularSystemAccessor msa(ms);

    for(int iMol=0;iMol<ms.MoleculesCount();iMol++){
        int nAtomsInMol = ms[iMol].AtomsCount();
        if(nAtomsInMol==0)
            continue;
        vector<int> globalIndexes;
        for(int iAtom=0;iAtom<nAtomsInMol;iAtom++){
            globalIndexes.push_back(msa.GlobalIndexOfAtom(ms[iMol][iAtom].globalSerial));
        }
        shuffle(globalIndexes.begin(),globalIndexes.end(),e);
        //Add a random chain to connect all atoms:
        for(int i=0;i<globalIndexes.size()-1;i++){
            int from  = globalIndexes[i];
            int to    = globalIndexes[i+1];
            if(not msa.IfBonded(from,to)){
                string fromSerial = msa.AtomByGlobalIndex(from).serial;
                string toSerial = msa.AtomByGlobalIndex(to).serial;
                Bond b = {fromSerial,toSerial,"pseudo",0};
                ms[iMol].bonds.push_back(make_shared<Bond>(b));
            }
        }
    }
}


// Test Randomly splitting a molecule, combine them, and find bonds
TEST_F(TestSplitting, split_and_find_bonds){

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
            cout<<"To "<<n<<" Piece(s):"<<endl;
            cout << "Splitted: \t"<<thems.Summary() << endl;
            // Then auto detect
            EXPECT_EQ(thems.BondsCount(),nBonds);
            cout << "AutoDete: \t"<<thems.Summary() << endl;

            // Random Translate
            thems.Translate(XYZ(urd(e),urd(e),urd(e)));
            thems.DetectBonds(GetDefaultBondDetector(5.0),true);
            EXPECT_EQ(thems.BondsCount(),nBonds);
            cout << "Shifted: \t"<<thems.Summary() << endl;

        }

        thems.Write(&mol2file, DATAFILESPATH + "/../dump" + to_string(iNo)+".mol2");
    }
}


// Test randomly split a molecule, delete all inter-molecular bonds, split the molsys into multiple molecules
TEST_F(TestSplitting, divide_by_connectivity){
    uniform_real_distribution<double> urd(-5,5);
    for(int iNo=1;iNo<ms.size();iNo++) {

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
            cout<<"Testing "<<n<<" pieces"<<endl;
            auto brokenMol = thems.DeepCopy();
            MolSysRandomSplit(brokenMol, n);
            cout<<"1. Random Split: "<<brokenMol.Summary()<<endl;
            // Delete inter-molecular bonds:
            brokenMol.interMolecularBonds.clear();
            // Artificially connect all atoms in each split molecules(no physical meaning)
            ConnectEveryAtomInEachMolecule(brokenMol);
            brokenMol.RenumberAtomSerials();
            // Find how many non-empty molecules
            int nMols = 0;
            int nBonds = brokenMol.BondsCount();
            int nAtoms = brokenMol.AtomsCount();
            for(int iMol=0;iMol<brokenMol.MoleculesCount();iMol++){
                if(brokenMol[iMol].AtomsCount()>0)
                    nMols++;
            }
            cout<<"2. Add Connect : "<<brokenMol.Summary()<<" : "<<nMols<<" non-empty mols"<<endl;
            // Then combine this artificial mol into a whole:
            MolSysReduceToSingleMolecule(brokenMol);
            // Then ask the algorithm to split the system according to connectivity, and expect the
            // correct number of molecules;
            MolSysSplitByConnectivity(brokenMol);
            cout<<"3. By connct   : "<<brokenMol.Summary()<<endl;
            EXPECT_EQ(brokenMol.MoleculesCount(), nMols);
            EXPECT_EQ(brokenMol.BondsCount(), nBonds);
            EXPECT_EQ(brokenMol.AtomsCount(),nAtoms);
        }

        // thems.Write(&mol2file, DATAFILESPATH + "/../dump" + to_string(iNo)+".mol2");
    }
}
