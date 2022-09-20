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

TEST(multithread,DISABLED_t1){
    while(true)
        TestMultiThread();
}

TEST(Reading,t1){
    MolecularSystem ms;
    QuickOpen(ms,curPath+"water.data");
    cout<<ms.Summary()<<endl;
    Trajectory traj(ms);
//    set<int> certainFrames;
//    for(int i=4;i<250;i+=3){
//        certainFrames.insert(i*10000);
//    }
    traj.Read(curPath+"water.lammpstrj");
}

TEST(KeywordsColumnPos,DISABLED_t1){
    KeywordsColumnPos kcp;
    kcp.FindColumnPos("ITEM: ATOMS id mol type x y z vx vy vz");
}