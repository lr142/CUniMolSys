#ifndef CUNIMOLSYS_VMDINTERFACE_H
#define CUNIMOLSYS_VMDINTERFACE_H
#include "utility.h"
#include "universalmolecularsystem.h"
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;

class VMDInterface{
    string molecularFileName;
    string scriptFileName;
    bool autobonds;
    double textsize;
    double textthickness;
    vector<string> commands;
    void Step1_Representations();
    void Step2_LoadMolecule();
    void Step3_Representations();
public:
    VMDInterface(string molecularFileName,string scriptFileName="vmd.tcl",bool autobonds=false);
    void AddLabels(vector<int> indexes,vector<string> labels,double *offset);
    void AddCommand(string cmd);
    void AddPeriodicBoundaryBox(Boundary boundary,string color="black",int width=1,string style="solid");
    void Run();
};
#endif //CUNIMOLSYS_VMDINTERFACE_H
