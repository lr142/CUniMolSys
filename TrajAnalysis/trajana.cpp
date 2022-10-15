#include "trajana.h"
#include "bonddetector.h"
#include "molecularmanipulator.h"
#include "lammpsdatafile.h"
#include "opener.h"
#include <algorithm>
using namespace std;


Analyzer::Analyzer(string path,int max_workers):max_workers(max_workers){
    Read(path);
}

void Analyzer::Read(string path){
    this->working_dir = fs::path(path);
    fs::path datafile = working_dir / "system.data";
    if(not fs::is_regular_file(datafile))
        ERROR("File "+datafile.string()+" doesn't exist!");
    this->ms = make_shared<MolecularSystem>("System");
    QuickOpen(*ms, datafile.string());
    // Create the accessor
    this->msa = make_shared<MolecularSystemAccessor>(*ms);
    cout<<ms->Summary()<<endl;

    // Create the trajectory, Find all lammpstrj files
    this->traj = make_shared<Trajectory>(*ms);
    auto itr = fs::directory_iterator(working_dir);
    vector<string> lammpstrajs;
    for(auto &item:itr){
        if(StringRegexMatch(item.path().filename().string(),"lammpstrj",false)){
            auto a = item.path().filename().string();
            lammpstrajs.push_back(a);
        }
    }
    // Sort before reading to assure correct order.
    sort(lammpstrajs.begin(),lammpstrajs.end());
    // Now read all trajectories
    for(auto &item:lammpstrajs){
        auto full_path = working_dir / item;
        traj->Read(full_path.string(),this->max_workers,100);
    }

    // output debugging info, optional
//    cout<<ms->AtomsCount()<<endl;
//    for(int i=0;i<traj->NFrames();i++) {
//        cout << "iFrame = " << i << " " << (*traj)[0].NAtoms() << endl;
//    }
//    for(int i=1;i<(*traj)[0].NAtoms();i++){
//        int diff = i - (*traj)[0].SerialToIndex(i);
//        ERROR_IF_FALSE(diff==1)
//    }
}

void Analyzer::locate_clay(vector<vector<int>>& atoms_indexes){
    for(int iMol=0;iMol<ms->MoleculesCount();iMol++){
        if((*ms)[iMol].AtomsCount()<8000)
            continue;
        vector<int> layer((*ms)[iMol].AtomsCount());
        for(int iAtom=0;iAtom<(*ms)[iMol].AtomsCount();iAtom++){
            layer[iAtom] = msa->GlobalIndexOfAtom((*ms)[iMol][iAtom]);
        }
        atoms_indexes.push_back(layer);
    }
}

void Analyzer::locate_polymers(vector<vector<int>>& atoms_indexes){
    set<string> acceptableElements = {"C","H","O","N","S"};
    vector<bool> mask(ms->AtomsCount(),false);
    // Polymers may be in one molecule, may be not.
    // Let work on a copy of ms, not the ms itself
    MolecularSystem ms_copy = ms->DeepCopy();
    // Let's at first find all polymer atoms, then extract all of them, and split into multiple chains.
    for(int iMol=0;iMol<ms_copy.MoleculesCount();iMol++) {
        int atomCount = ms_copy[iMol].AtomsCount();
        if (atomCount >= 8000 or atomCount <= 3)
            continue;
        bool isPolymer = true;
        for (int iAtom = 0; iAtom < atomCount; iAtom++) {
            auto &atom = ms_copy[iMol][iAtom];
            if (not acceptableElements.contains(atom.element)) {
                isPolymer = false;
                break;
            }
        }
        if (not isPolymer)
            continue;
        for (int iAtom = 0; iAtom < atomCount; iAtom++) {
            auto &atom = ms_copy[iMol][iAtom];
            int index = msa->GlobalIndexOfAtom(atom);
            // before we process, let's at first record the atoms's original index in the atom, in the parent property
            atom.parent = to_string(index);
            mask[index] = true;
        }
    }
    // Then select out the polymer atoms, and split them into multiple chains
    MolSysSubsystemByMask(ms_copy,mask);
    ms_copy.DetectBonds(GetDefaultBondDetector(5.0),true);
    MolSysSplitByConnectivity(ms_copy);
    // Record the results
    for(int iMol=0;iMol<ms_copy.MoleculesCount();iMol++){
        vector<int> indexes(ms_copy[iMol].AtomsCount());
        for(int iAtom=0;iAtom<ms_copy[iMol].AtomsCount();iAtom++){
            // the index is previously recorded in 'parent' property
            indexes[iAtom] = stoi(ms_copy[iMol][iAtom].parent);
        }
        sort(indexes.begin(),indexes.end());
        atoms_indexes.push_back(indexes);
    }
}
void Analyzer::locate_ions(vector<vector<int>>& atoms_indexes, std::set<string> &ions){
    vector<int> indexes;
    for(int iMol=0;iMol<ms->MoleculesCount();iMol++){
        if( (*ms)[iMol].AtomsCount() >= 8000 )
            continue;
        for(int iAtom=0;iAtom<(*ms)[iMol].AtomsCount();iAtom++){
            auto &atom = (*ms)[iMol][iAtom];
            string element = atom.element;
            if( not ions.contains(element))
                continue;
            indexes.push_back(msa->GlobalIndexOfAtom(atom));
        }
    }
    atoms_indexes.push_back(indexes);
}
void Analyzer::locate_anions(vector<vector<int>>& atoms_indexes){
    set<string> anions = {"Cl"};
    locate_ions(atoms_indexes,anions);
}
void Analyzer::locate_cations(vector<vector<int>>& atoms_indexes){
    set<string> cations = {"K","Na","Ca"};
    locate_ions(atoms_indexes,cations);
}
void Analyzer::locate_water(vector<vector<int>>& atoms_indexes){
    set<string> acceptableElements = {"O","H"};
    for(int iMol=0;iMol<ms->MoleculesCount();iMol++){
        if((*ms)[iMol].AtomsCount()!=3)
            continue;
        bool waterFlag = true;
        vector<int> this_water(3);
        for(int iAtom=0;iAtom<3;iAtom++){
            auto &atom = (*ms)[iMol][iAtom];
            this_water[iAtom] = msa->GlobalIndexOfAtom(atom);
            if(not acceptableElements.contains(atom.element))
                waterFlag = false;
        }
        if(waterFlag)
            atoms_indexes.push_back(this_water);
    }
}
void Analyzer::Locate(vector<vector<int>> &atoms_indexes, ComponentTypes comp) {
    switch(comp){
        case POLYMER:
            locate_polymers(atoms_indexes);break;
        case CLAY:
            locate_clay(atoms_indexes);break;
        case ANIONS:
            locate_anions(atoms_indexes);break;
        case CATIONS:
            locate_cations(atoms_indexes);break;
        case WATER:
            locate_water(atoms_indexes);break;
        default:
            ERROR("Invalid ComponentType");
    }
}

fs::path Analyzer::get_output_dir() {
    fs::path output_dir = working_dir / "analyze";
    if (not fs::exists(output_dir))
        fs::create_directories(output_dir);
    else if(not fs::is_directory(output_dir)){
        fs::remove_all(output_dir);
        fs::create_directories(output_dir);
    }
    return output_dir;
}

double Analyzer::calc_clay_top() {
    vector<vector<int>> clays;
    Locate(clays,CLAY);
    double clayTop = -MY_LARGE;
    for(int i=0;i<clays[0].size();i++){
        int index = clays[0][i];
        double z = msa->AtomByGlobalIndex(index).xyz[2];
        clayTop = max(clayTop,z);
    }
    return clayTop;
}

void Analyzer::PolymerZPositions() {

    double clayTop = calc_clay_top();
    vector<vector<int>> polymers;
    Locate(polymers,POLYMER);

    // Calculate positions;
    int NPolymers = (int)polymers.size();
    int NFrames = traj->NFrames();
    vector<vector<double>> pos_avg(NFrames);
    vector<vector<double>> pos_min(NFrames);

    for(int iFrame=0;iFrame<NFrames;iFrame++){
        pos_avg[iFrame].resize(NPolymers);
        pos_min[iFrame] = vector<double>(NPolymers,MY_LARGE);
        for(int iPoly=0;iPoly<NPolymers;iPoly++){
            int NAtoms = (int)polymers[iPoly].size();
            double sum_of_z = 0.0;
            for(int iAtom=0;iAtom<NAtoms;iAtom++){
                int index = polymers[iPoly][iAtom];
                double z = (*traj)[iFrame].X()[index][2];
                sum_of_z += z;
                pos_min[iFrame][iPoly] = min(pos_min[iFrame][iPoly],z );
            }
            double avg_z = sum_of_z / NAtoms;
            avg_z -= clayTop;
            pos_avg[iFrame][iPoly] = avg_z;
            pos_min[iFrame][iPoly] -= clayTop;
        }
    }


    // output
    ofstream ofs[2];
    ofs[0].open(get_output_dir() / "PolymerZ_avg.csv");
    ofs[1].open(get_output_dir() / "PolymerZ_min.csv");
    for(int iFile=0; iFile<2; iFile++) {
        ofs[iFile] << "#, NRows (one line for each frame); NColumns (1st column for iFrame, then one column for each chain)"
            << endl;
        ofs[iFile] << NFrames << "," << NPolymers + 1 << endl;
        for (int i = 0; i < NFrames; i++) {
            ofs[iFile] << i << ",";
            for (int j = 0; j < NPolymers; j++) {
                ofs[iFile] << (iFile==0 ? pos_avg[i][j] : pos_min[i][j]) << ",";
            }
            ofs[iFile] << endl;
        }
    }
}

int main(int argc,char* argv[]){
    string path;
    int max_workers;
    if(argc<2)
        path = "./";
    else
        path = string(argv[1]);
    max_workers = argc<3 ? -1 : atoi(argv[2]);
    Analyzer ana(path,max_workers);

    ana.PolymerZPositions();
}
