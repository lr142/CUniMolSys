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

enum ComponentTypes {CLAY,POLYMER,ANIONS,CATIONS,WATER,COO,AMMONIUM,CONH2};

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

    /* Locate a functional group from a set of atoms defined in molecule.
    * write the results in 'found'. size() of 'found' is the number of functional groups found.
    * centerElement is the element at the center of the functional group
     * bondCount is the number of bonds attached to the center atom
     * bondedElement is the elements that bonded to the center atom and considered as part of the functions group.
     * size of 'bondedElement' must be <= 'bondCount' */
    void locate_functional_group(vector<int> &molecule, vector<vector<int>> &found,
                                 string centerElement, int bondCount, vector<string> bondedElements);
    void locate_COO(vector<vector<int>>& atoms_indexes);
    void locate_AmmoniumN(vector<vector<int>>& atoms_indexes);
    void locate_CONH2(vector<vector<int>>& atoms_indexes);

    fs::path get_output_dir();
    double calc_clay_top();

    vector<double> rdf_between(vector<vector<int>> &from, vector<vector<int>> &to, double range_in_nm, double binSize_in_nm);
    vector<double> msd_of(vector<int> &atoms);
public:
    Analyzer(string path="./",int max_workers=-1);
    ~Analyzer() = default;
    void Read(string path);
    /* See the definition of locate_xxx above */
    void Locate(vector<vector<int>>& atoms_indexes, ComponentTypes comp);
    // Draw out the selected subsystem.
    void Locate_Demonstrate(vector<vector<int>>& atoms_indexes,string filename);
    void PolymerZPositions();
    void ParticleZDistribution(double binSize_in_nm=0.02);
    void RDF(double range_in_nm=3.0, double binSize_in_nm=0.02);
    void MSD();
    void Density();
    void ChainLengths(double range_in_nm=3.0, double binSize_in_nm=0.02);
};

void MainAsAnalyzePolymerDesorptionZ(int argc,char* argv[]);

// Unit Tests:
void __TestLocateModule__(int argc,char* argv[]);


#endif //CUNIMOLSYS_TRAJANA_H
