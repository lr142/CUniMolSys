#include <iostream>
#include "gtest/gtest.h"
#include "opener.h"
#include "molecularmanipulator.h"
#include "universalmolecularsystem.h"
#include "trajectory.h"
#include <vector>
#include <fstream>

using namespace std;

void TestMultiThread(){
    MolecularSystem ms;
//    QuickOpen(ms,DATAFILESPATH+"/Structures/tip4p_water_5nm_chunk.mol2");
    Trajectory traj(ms);
//    cout<<"NAtoms = "<<ms.AtomsCount()<<endl;
    traj._testMultiThread();
}

string curPath = DATAFILESPATH+"/../UnitTests/008.TestTraj/";

TEST(KeywordsColumnPos,DISABLED_t1){
    KeywordsColumnPos kcp;
    kcp.FindColumnPos("ITEM: ATOMS id mol type x y z vx vy vz");
}

TEST(multithread,DISABLED_t2){
    for(;;) {
        TestMultiThread();
    }
}

TEST(Reading,DISABLED_t1){
    MolecularSystem ms;
    QuickOpen(ms,curPath+"water.data");
    cout<<ms.Summary()<<endl;
    Trajectory traj(ms);
//    set<int> certainFrames;
//    for(int i=4;i<250;i+=3){
//        certainFrames.insert(i*10000);
//    }
    traj.Read(curPath+"water.traj");
}

TEST(Reading,DISABLED_t2){
    MolecularSystem ms;
    curPath = DATAFILESPATH+"/../UnitTests/008.TestTraj/";
    QuickOpen(ms,curPath+"polymer.data");
    cout<<ms.Summary()<<endl;
    if(false){
        MolecularSystem copy = ms.DeepCopy();
        MolSysReduceToSingleMolecule(copy);
        QuickSave(copy, DATAFILESPATH + "/../dump.mol2");
        exit(0);
    }
    {
        Trajectory traj(ms);
//            set<int> certainFrames;
//            for(int i=4;i<250;i+=3){
//                certainFrames.insert(i*10000);
//            }
        for(int i=0;i<10;i++) {
            traj.Read(curPath + "polymer.traj");
            //traj.ShowTrajectory(DATAFILESPATH+"/../dump.mol2",false);
        }
    }
}


TEST(Reading,largeSys){
    MolecularSystem ms;
    curPath = DATAFILESPATH+"/../UnitTests/008.TestTraj/";
    QuickOpen(ms,curPath+"polymer.data");
    cout<<ms.Summary()<<endl;
    if(false){
        MolecularSystem copy = ms.DeepCopy();
        MolSysReduceToSingleMolecule(copy);
        QuickSave(copy, DATAFILESPATH + "/../dump.mol2");
        exit(0);
    }
    {
        Trajectory traj(ms);
        set<int> certainFrames;
//        for(int i=0;i<6000000;i+=100000){
//            certainFrames.insert(i);
//        }
        traj.Read(curPath + "system.lammpstrj",-1,MY_LARGE,true,certainFrames);
        traj.Read(curPath + "system.lammpstrj.2",-1,MY_LARGE,true,certainFrames);
        for(int i=0;i<traj.NFrames();i++){
            cout<<"Frame = "<<i<<", ts = "<<traj[i].ts_<<", NAtoms = "<<traj[i].nAtoms_<<endl;
        }
//        traj.ShowTrajectory(DATAFILESPATH+"/../dump.mol2",false);
    }
}

