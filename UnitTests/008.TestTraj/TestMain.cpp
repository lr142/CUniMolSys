#include <iostream>
#include "gtest/gtest.h"
#include "opener.h"
#include "molecularmanipulator.h"
#include "universalmolecularsystem.h"
#include "trajectory.h"
#include <vector>
#include <fstream>

using namespace std;

string curPath = DATAFILESPATH+"/../UnitTests/008.TestTraj/";

TEST(KeywordsColumnPos,DISABLED_t1){
    KeywordsColumnPos kcp;
    kcp.FindColumnPos("ITEM: ATOMS id mol type x y z vx vy vz");
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
    QuickSave(ms,DATAFILESPATH+"/../oridump.mol2");
    if(false){
        MolecularSystem copy = ms.DeepCopy();
        MolSysReduceToSingleMolecule(copy);
        QuickSave(copy, DATAFILESPATH + "/../dump.mol2");
        exit(0);
    }
    {
        Trajectory traj(ms);
        set<long long> certainFrames;
//        for(int i=4;i<100;i+=2){ // 4,6,8,10,12,... * 1000
//            certainFrames.insert(i*1000);
//        }
        for(int i=0;i<1;i++) {
            traj.Read(curPath + "polymer.traj",-1,99999,true,certainFrames); // 4,6,8,10,12 * 10000
            for(int i=0;i<traj.NFrames();i++){
                cout<<"Frame = "<<i<<", ts = "<<traj[i].TS()<<", NAtoms = "<<traj[i].NAtoms()<<endl;
            }
            set<int> update_atoms;
            for(int i=16500;i<18000;i++){
                update_atoms.insert(i);
            }
            set<int> include_atoms;
            for(int i=16500;i<19000;i++){
                include_atoms.insert(i);
            }
            traj.ShowTrajectory(DATAFILESPATH+"/../dump.mol2",update_atoms,include_atoms,4);
        }
    }
}


TEST(Reading,largeSys){
    MolecularSystem ms;
    curPath = DATAFILESPATH+"/../../0402_PolyDADMAC_400K/";
    QuickOpen(ms,curPath+"system.data");
    cout<<ms.Summary()<<endl;
    {
        Trajectory traj(ms);
        set<long long> certainFrames;
        for(int i=0;i<500;i+=4){
            certainFrames.insert(i*10000);
        }
        traj.Read(curPath + "system.lammpstrj",4,99999,true);//,certainFrames);
        traj.Read(curPath + "system.lammpstrj.2",4,99999,true);
//        for(int i=0;i<traj.NFrames();i++){
//            cout<<"Frame = "<<i<<", ts = "<<traj[i].TS()<<", NAtoms = "<<traj[i].NAtoms()<<endl;
//        }
        string filename = DATAFILESPATH+"/../dump.mol2";
        set<int> update;
        for(int i=16500;i<18000;i++)
            update.insert(i);
        traj.ShowTrajectory(filename,update,update,4);
    }
}

