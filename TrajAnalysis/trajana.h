#ifndef CUNIMOLSYS_TRAJANA_H
#define CUNIMOLSYS_TRAJANA_H

#include "universalmolecularsystem.h"
#include "trajectory.h"
#include <iostream>
#include <string>
#include <memory>
#include <set>
#include <filesystem>

using std::vector;
using std::string;
namespace fs = std::filesystem;

enum ComponentTypes {CLAY,POLYMER,ANIONS,CATIONS,WATER};

class Analyzer{
    int max_workers;
    fs::path working_dir;
    shared_ptr<MolecularSystem> ms;
    shared_ptr<Trajectory> traj;
    shared_ptr<MolecularSystemAccessor> msa;
    /* A series of functions locating components is the molecule
     * These functions will add 1 or more vectors after the last element in atoms_indexes.
     * eg. locate_clay will append two lists, representing the lower and upper plates. locate_polymers will
     * append N lists, where N is the number of chains in the system.
     * In all lists, atoms index is the index of atom in the global system, starting from 0. */
    void locate_clay(vector<vector<int>>& atoms_indexes);
    void locate_polymers(vector<vector<int>>& atoms_indexes);
    void locate_ions(vector<vector<int>>& atoms_indexes, std::set<string> &ions);
    void locate_anions(vector<vector<int>>& atoms_indexes);
    void locate_cations(vector<vector<int>>& atoms_indexes);
    void locate_water(vector<vector<int>>& atoms_indexes);
    fs::path get_output_dir();
    double calc_clay_top();

public:
    Analyzer(string path="./",int max_workers=-1);
    ~Analyzer() = default;
    void Read(string path);
    /* See the definition of locate_xxx above */
    void Locate(vector<vector<int>>& atoms_indexes, ComponentTypes comp);
    void PolymerZPositions();
};

#endif //CUNIMOLSYS_TRAJANA_H
