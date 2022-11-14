#include "trajana.h"
#include "bonddetector.h"
#include "molecularmanipulator.h"
#include "lammpsdatafile.h"
#include "opener.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
using namespace std;

// Some common tool functions.
/* Get the geometry center of atoms in frame iFrame, atom indexes are in the vector atoms.
 * Must conform the PBC defined in U,V,W. This function does not rely on the I() 'image flags' info in trajectory */
XYZ __centerPosition__(Trajectory &traj, vector<int>& atoms, int iFrame, double U,double V,double W){
    auto xyzs = traj[iFrame].X();
    XYZ UVW(U,V,W);
    XYZ center;
    XYZ firstAtom = xyzs[atoms[0]];
    center = firstAtom;
    if(atoms.size()==1)
        return center;
    // Must get the unwrapped coordinates, by checking the distance to the 1st atom, must < L/2
    for(int i=1;i<atoms.size();i++){
        XYZ curAtom = xyzs[atoms[i]];
        for(int j=0;j<3;j++){
            while(curAtom[j]-firstAtom[j] > UVW[j]/2.0)
                curAtom[j] -= UVW[j];
            while(firstAtom[j] - curAtom[j] > UVW[j]/2.0)
                curAtom[j] += UVW[j];
        }
        center += curAtom;
    }
    center *= (1.0/atoms.size());
    return center;
}
/* Distance of two atoms in PBC. It does not rely on the image flags I(), relies on PBC info instead */
double __dist_in_PBC__(XYZ& p1, XYZ& p2, double U,double V, double W){
    double bound[3] = {U,V,W};
    double squared_sum = 0.0;
    for(int i=0;i<3;i++){
        double diff = fabs(p1[i] - p2[i]);
        while(diff > bound[i])
            diff -= bound[i];
        diff = min(diff, bound[i]-diff);
        squared_sum += pow(diff, 2);
    }
    return sqrt(squared_sum);
}
/* "flatten" a vector of vectors into a 1d vector */
vector<int> __flatten__(vector<vector<int>>& original){
    vector<int> result;
    for(auto &v:original){
        for(auto &i:v){
            result.push_back(i);
        }
    }
    return result;
}
/* Get the unwrapped coordinate of atom iAtom in a certain trajectory frame */
XYZ __unwrapped_XYZ_(TrajectoryFrame& frame,int iAtom,XYZ &UVW){
    XYZ xyz = frame.X()[iAtom];
    for(int i=0;i<3;i++){
        xyz[i] += frame.I()[iAtom][i] * UVW[i];
    }
    return xyz;
}
/* A 1-d linear regression. Write the result into slope and intercept. */
void __linear_fitting__(vector<double> &data, int start, int end, double &slope, double &intercept){
    double sx = 0.0,sy = 0.0,sxx = 0.0,sxy = 0.0;
    int n = end-start;
    for(int x=start;x<end;x++){
        double y = data[x];
        sx += x;
        sy += y;
        sxx += x*x;
        sxy += x*y;
    }
    intercept = (sy*sxx - sx*sxy) / ( n*sxx - sx*sx);
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx);
}

Analyzer::Analyzer(string path,int max_workers):max_workers(max_workers){
    Read(path);
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
        traj->Read(full_path.string(),this->max_workers);
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
        case COO:
            locate_COO(atoms_indexes);break;
        case CONH2:
            locate_CONH2(atoms_indexes);break;
        case AMMONIUM:
            locate_AmmoniumN(atoms_indexes);break;
        default:
            ERROR("Invalid ComponentType");
    }
}
void Analyzer::locate_functional_group(vector<int> &molecule, vector<vector<int>> &found,
                                       string centerElement, int bondCount, vector<string> bondedElements){
    for(int iAtom:molecule){
        auto &centerAtom = msa->AtomByGlobalIndex(iAtom);
        vector<int> bondedList;
        for(auto &item:msa->GetBondedMap()[iAtom]){
            // item.first is the bonded atom index,
            // item.second points to the actual Bond structure.
            bondedList.push_back(item.first);
        }
        if(centerAtom.element!=centerElement or bondedList.size()!=bondCount)
            continue;
        vector<int> group;
        group.push_back(iAtom);
        // check all required bonded elements are present. In case some elements will appear more than once, set a flag
        // to note that this atom is already added
        vector<bool> visited(bondCount,false);
        bool groupMatched = true;  // this is set to false if any atom is not matched.
        for(string expected_element:bondedElements){
            bool atomMatched = false; // this is set to true if any candidate atom can match this atom
            for(int i=0;i<bondCount;i++){
                if(visited[i])
                    continue;
                int iCandidateAtom = bondedList[i];
                if(msa->AtomByGlobalIndex(iCandidateAtom).element == expected_element){
                    visited[i] = atomMatched = true;
                    group.push_back(iCandidateAtom);
                    break;
                }
            }
            if(not atomMatched) {
                groupMatched = false;
                break;
            }
        }
        if(groupMatched){
            found.push_back(group);
        }
    }
}
void Analyzer::locate_COO(vector<vector<int>>& atoms_indexes){
    vector<vector<int>> list;
    locate_polymers(list);
    for(auto &mol: list) {
        locate_functional_group(mol, atoms_indexes, "C", 3, vector<string>{"O", "O"});
    }
}
void Analyzer::locate_AmmoniumN(vector<vector<int>>& atoms_indexes){
    vector<vector<int>> list;
    locate_polymers(list);
    for(auto &mol: list) {
        locate_functional_group(mol, atoms_indexes, "N", 4, vector<string>{});
    }
}
void Analyzer::locate_CONH2(vector<vector<int>>& atoms_indexes){
    vector<vector<int>> list;
    locate_polymers(list);
    vector<vector<int>> CO;
    vector<vector<int>> NH2;
    for(auto &mol: list) {
        locate_functional_group(mol, CO, "C", 3, vector<string>{"O","N"});
        locate_functional_group(mol, NH2, "N", 3, vector<string>{"H","H"});
    }
    // now we must combine the two lists: CO and NH2. since C is bonded to N, we can search the last item of each list in CO:
    for(auto &co_group: CO){
        int Nindex = co_group[2]; // 0,1,2 are C,O,and N, respectively
        for(auto &nh2_group:NH2){
            if(Nindex == nh2_group[0]){
                // found a match, these two are the same group. Totally 5 atoms
                atoms_indexes.push_back(vector<int>{co_group[0],co_group[1],co_group[2],nh2_group[1],nh2_group[2]});
                break;
            }
        }
    }
}
// Draw out the selected subsystem.
void Analyzer::Locate_Demonstrate(vector<vector<int>>& atoms_indexes,string filename){
    vector<bool> mask(msa->AtomsCount(),false);
    int atomCount = 0;
    for(auto &i:atoms_indexes){
        for(auto &j:i){
            mask[j] = true;
        }
        atomCount += i.size();
    }
    cout<<"Mols: "<<atoms_indexes.size()<<", Atoms: "<<atomCount<<endl;

    auto ms_copy = ms->DeepCopy();
    MolSysSubsystemByMask(ms_copy,mask);
    QuickSave(ms_copy,filename);
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

void Analyzer::ParticleZDistribution(double binSize_in_nm){
    // find what are in the system
    map<string, vector<vector<int>>> components;
    vector<vector<int>> temp;
    //temp.clear();locate_clay(temp);components["CLAY"] = temp;
    //temp.clear();locate_polymers(temp);components["POLYMER"] = temp;
    temp.clear();locate_anions(temp);components["ANIONS"] = temp;
    temp.clear();locate_cations(temp);components["CATIONS"] = temp;
    temp.clear();locate_AmmoniumN(temp);components["AMMONIUM"] = temp;
    temp.clear();locate_COO(temp);components["COO"] = temp;
    temp.clear();locate_CONH2(temp);components["CONH2"] = temp;
    double clay_top = this->calc_clay_top();
    // for Anions and Cations, originally there is only 1 vector in the 1st dimension, now we split it into nIons vecs,
    // each element of which is a vector of length 1
    for(string key:vector<string>{"ANIONS","CATIONS"}){
        int nIons= components[key][0].size();
        vector<vector<int>> new_vec(nIons);
        for(int i=0;i<nIons;i++){
            new_vec[i] = vector<int>{components[key][0][i]};
        }
        components[key] = new_vec;
    }
    // Now we create a histogram for each type of particles, count the distribution of each type of particles in each bin
    // Current set the range to be 0nm to 4 nm above the clay
    int NBins = (int)ceil(4.0 / binSize_in_nm);
    int NFrames = this->traj->NFrames();
    vector<string> keys = {"CATIONS","ANIONS","COO","AMMONIUM","CONH2"};
    map<string,vector<double>> histograms;
    for(auto key:keys){
        histograms[key] = vector<double>(NBins,0.0);
        vector<vector<int>> &groupsIndexes = components[key];
        for(int iFrame=0;iFrame<NFrames;iFrame++) {
            auto posMatrix = (*traj)[iFrame].X();
            for (vector<int> &group: groupsIndexes) {
                // Get the average coord of this group
                int NAtomsInThisGroup = group.size(); // should be 1 for cations/ions/Ammonium. 3 for COO and 5 for CONH2
                XYZ avg_coord = XYZ(0, 0, 0);
                for (int iAtom: group) {
                    avg_coord += posMatrix[iAtom];
                }
                avg_coord *= 1.0 / NAtomsInThisGroup;

                // Check which bin it belongs by the Z-coord
                double dist_in_nm = (avg_coord[2] - clay_top) / 10.0;
                int binIndex = (int) floor(dist_in_nm / binSize_in_nm);
                if (binIndex < 0 or binIndex >= NBins) // out of range
                    continue;
                histograms[key][binIndex] += 1.0;
            }
        }
        // Calculate frames average
        for(auto &val:histograms[key])
            val /= NFrames;
        // Convert counts into densities. Density unit in counts/nm^3, then additional to mol/L
        double vol = ms->boundary.GetU()[0] * ms->boundary.GetV()[1] * 0.01 * binSize_in_nm;
        for(auto &val:histograms[key]) {
            val /= vol;   // unit in count/nm^3
            val *= (10/6.02214);   // 1E24 / 6.02214E23;      // unit in mol/L
        }
    }
    // output
    ofstream ofs;
    ofs.open(get_output_dir() / "DistributionZ.csv");
    ofs << "#, NRows (dist), NColumns"<<endl;
    ofs << "#, ";
    for(auto key:keys)
        ofs<<key<<", ";
    ofs<<endl;
    ofs<<NBins<<","<<keys.size()+1<<endl;
    for (int i = 0; i < NBins; i++) {
        ofs << (i+0.5)*binSize_in_nm << ",";
        for (int j = 0; j < keys.size(); j++) {
            ofs<<histograms[keys[j]][i]<<",";
        }
        ofs << endl;
    }
}

vector<double> Analyzer::rdf_between(vector<vector<int>> &from, vector<vector<int>> &to, double range_in_nm, double binSize_in_nm){
    int nAtomsA = from.size();
    int nAtomsB = to.size();
    int nBins = (int)ceil(range_in_nm/binSize_in_nm);

    double U = ms->boundary.GetU()[0];
    double V = ms->boundary.GetV()[1];
    double W = ms->boundary.GetW()[2];

    vector<double> counter(nBins,0.0);
    // counter the neighbors of A (which is B) for each A in each frame.

    // for debuggin, check clashing atoms
    //vector<vector<string>> dump;

    for(int iFrame=0;iFrame<traj->NFrames();iFrame++) {

//        dump.push_back(vector<string>());

        for (int i = 0; i < nAtomsA; i++) {
            for (int j = 0; j < nAtomsB; j++) {
                XYZ xyzFrom = __centerPosition__(*traj, from[i], iFrame,U,V,W);
                XYZ xyzTo = __centerPosition__(*traj, to[j], iFrame,U,V,W);
                double dist = __dist_in_PBC__(xyzFrom, xyzTo, U, V, W) * 0.1;  // also convert distance to nm
//                if(dist < 0.1) {
//                    for(auto iFrom:from[i]){
//                        ostringstream oss;
//                        oss<<msa->AtomByGlobalIndex(iFrom).element<<" ";
//                        auto &xyz = (*traj)[iFrame].X()[iFrom];
//                        oss<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2];
//                        dump[dump.size()-1].push_back(oss.str());
//                    }
//                    for(auto iTo:to[j]){
//                        ostringstream oss;
//                        oss<<msa->AtomByGlobalIndex(iTo).element<<" ";
//                        auto &xyz = (*traj)[iFrame].X()[iTo];
//                        oss<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2];
//                        dump[dump.size()-1].push_back(oss.str());
//                    }
//                }

                int whichBin = (int) floor(dist / binSize_in_nm);
                ERROR_IF_TRUE(whichBin < 0)
                if (whichBin < nBins)
                    counter[whichBin] += 1.0;
            }
        }
    }

    // debugging for checking clash
//    ofstream dumpfile("dump.xyz");
//    for(auto &frame:dump){
//        if(frame.size()==0)
//            continue;
//        else
//            dumpfile<<frame.size()<<endl<<endl;
//        for(auto &line:frame){
//            dumpfile<<line<<endl;
//        }
//    }

    // calculate the average neighbors per group A per frame
    for(int i=0;i<nBins;i++){
        counter[i] /= (nAtomsA * traj->NFrames());
    }
    // calculate the un-normalized RDF by divide the numbers by volumes of the shell
    vector<double> rdf_un_ormalized(nBins);
    double total_count = 0;
    for(int i=0;i<nBins;i++){
        double radius = (i+0.5)*binSize_in_nm;
        double vol_of_shell = 4 * MY_PI * pow(radius,2) * binSize_in_nm;
        rdf_un_ormalized[i] = counter[i]/vol_of_shell;
        total_count += counter[i];
    }
    // calculate the normalized RDF by divide each number by the average density
    vector<double> rdf(nBins);
    double average_density = total_count / (4 * MY_PI * pow(range_in_nm,3) / 3);
    for(int i=0;i<nBins;i++){
        rdf[i] = rdf_un_ormalized[i]/average_density;
    }
    return rdf;
}
void Analyzer::RDF(double range_in_nm, double binSize_in_nm){
    map<string, vector<vector<int>>> components;
    vector<vector<int>> temp;
    //temp.clear();locate_clay(temp);components["CLAY"] = temp;
    //temp.clear();locate_polymers(temp);components["POLYMER"] = temp;
    temp.clear();locate_anions(temp);components["ANIONS"] = temp;
    temp.clear();locate_cations(temp);
    // for cations, exclude counter Na+ below the plates
    components["CATIONS"] = vector<vector<int>>(1);
    for(int i:temp[0]){
        if(msa->AtomByGlobalIndex(i).element=="Na")
            continue;
        else
            components["CATIONS"][0].push_back(i);
    }
    temp.clear();locate_AmmoniumN(temp);components["AMMONIUM"] = temp;
    temp.clear();locate_COO(temp);components["COO"] = temp;
    temp.clear();locate_CONH2(temp);components["CONH2"] = temp;

    // for Anions and Cations, originally there is only 1 vector in the 1st dimension, now we split it into nIons vecs,
    // each element of which is a vector of length 1
    for(string key:vector<string>{"ANIONS","CATIONS"}){
        int nIons= components[key][0].size();
        vector<vector<int>> new_vec(nIons);
        for(int i=0;i<nIons;i++){
            new_vec[i] = vector<int>{components[key][0][i]};
        }
        components[key] = new_vec;
    }
    vector<string> fromKeys = {"COO","AMMONIUM","CONH2"};
    vector<string> toKeys = {"CATIONS","ANIONS","AMMONIUM"};
    map<string,vector<double>> all_rdfs;
    for(auto fromKey:fromKeys){
        for(auto toKey:toKeys){
            auto &fromSet = components[fromKey];
            auto &toSet = components[toKey];
            if(fromSet.size()==0 or toSet.size()==0)
                continue;
            vector<double> rdf = rdf_between(fromSet,toSet,range_in_nm,binSize_in_nm);
            all_rdfs[fromKey+"--"+toKey] = rdf;
        }
    }
    // output
    ofstream ofs;
    ofs.open(get_output_dir() / "RDF.csv");
    ofs << "#, NRows (dist), NColumns, Details on the 2nd Line "<<endl;
    int nSeries = all_rdfs.size();
    int nBins = 0;
    ofs<<"#,";
    for(auto &item:all_rdfs){
        ofs<<item.first<<", ";
        nBins = item.second.size();
    }ofs<<endl;
    ofs<<nBins<<","<<nSeries+1<<endl;
    for (int i = 0; i < nBins; i++) {
        ofs<<(0.5+i)*binSize_in_nm<<", ";
        for(auto &item:all_rdfs){
            ofs<<item.second[i]<<", ";
        }
        ofs<<endl;
    }
}


vector<double> Analyzer::msd_of(vector<int> &atoms) {
    int nFrames = this->traj->NFrames();
    int nAtoms = atoms.size();
    vector<double> msd_values(nFrames-1, 0.0);
    // Calculate the MSD ... The process of calculating MSD is
    // SD(iAtom, iFrame_start, iFrame_end) = | r(iAtom,iFrame_start) - r(iAtom,iFrame_end)|^2
    // MSD(iFrame_start, iFrame_end) = Sigma_{all iAtoms}( SD(iAtom, iFrame_start, iFrame_end) )  / nAtoms;
    // MSD( Dframes) =  Sigma_{all iFrame_start, iFrame_end such that iFrame_start - iFrame_end == Dframes }( MSD(iFrame_start, iFrame_end) ) / N of such combinations
    // This means that the time to calculate MSD(short time interval) is proportional to the nFrames/DFrames.
    XYZ UVW = {ms->boundary.GetU()[0], ms->boundary.GetV()[1], ms->boundary.GetW()[2]};
    for(int dFrames=1; dFrames<nFrames-1; dFrames++){ // dFrames is the index of the result msd_values
        int nEvaluations = nFrames-dFrames; // nEvaluations is the number of ways to sample startFrame and endFrame such as endFrame-startFrame == dFrames;
        double MSD_evaluation_sum = 0.0;
        for(int iFrameStart=0; iFrameStart<nEvaluations; iFrameStart++ ){
            int iFrameEnd = iFrameStart+dFrames;
            double totalSD = 0.0;
            auto &frameEnd = traj->operator[](iFrameEnd);
            auto &frameStart = traj->operator[](iFrameStart);
            for(int i=0; i<nAtoms; i++){
                int iAtom = atoms[i];
                XYZ from = __unwrapped_XYZ_(frameStart,iAtom,UVW);
                XYZ to   = __unwrapped_XYZ_(frameEnd,iAtom,UVW);
                XYZ diff = from - to;
                double SD = diff.NormSquared();
                totalSD += SD;
            }
            totalSD /= nAtoms;
            MSD_evaluation_sum += totalSD;
        }
        msd_values[dFrames] =  MSD_evaluation_sum / nEvaluations;
    }
    // msd_values[0]用来记录斜率。只用 1/4到3/4这段数据来拟合。拟合得到的斜率与截距写入msd_values[]的最后两个元素
    double intercept, slope;
    __linear_fitting__(msd_values,nFrames/4,nFrames/4*3,slope,intercept);
    msd_values.push_back(slope);
    msd_values.push_back(intercept);
    // last element is self-diffusion coefficient (in unit of 10^-9 m^2/s)
    msd_values.push_back(slope/6.0);
    return msd_values;
}
void Analyzer::MSD(){
    map<string, vector<int>> components;
    vector<vector<int>> temp;
    // Locate the following components: ANIONS, CATIONS, POLYMER, WATER
    temp.clear();locate_anions(temp);components["ANIONS"] = __flatten__(temp);
    temp.clear();locate_cations(temp);
    // for cations, exclude counter Na+ below the plates
    components["CATIONS"] = vector<int>();
    for(int i:temp[0]){
        if(msa->AtomByGlobalIndex(i).element=="Na")
            continue;
        else
            components["CATIONS"].push_back(i);
    }
    temp.clear();locate_polymers(temp);components["POLYMER"] = __flatten__(temp);
    temp.clear();locate_water(temp);components["WATER"] = __flatten__(temp);

    map<string,vector<double>> MSD_values;

    vector<string> keys;
    for(auto &item:components){
        auto key = item.first;
        auto atoms = item.second;
        MSD_values[key] = msd_of(atoms);
        keys.push_back(key);
    }

    // output
    ofstream ofs;
    ofs.open(get_output_dir() / "MSD.csv");
    ofs << "#, NRows (dist), NColumns, Details on the 2nd Line "<<endl;
    int nSeries = keys.size();
    int nLines = MSD_values[keys[0]].size();
    ofs<<"#,";
    for(auto &k:keys)
        ofs<<k<<",";
    ofs<<endl<<nLines<<","<<nSeries+1<<endl;
    for (int i = 0; i < nLines; i++) {
        ofs<<i<<", ";
        for(auto &k:keys){
            ofs<<MSD_values[k][i]<<", ";
        }
        ofs<<endl;
    }
}
void Analyzer::Density(){
    vector<vector<int>> clay;
    this->locate_clay(clay);
    double clay_lower = calc_clay_top();
    int clay_upper_index = -1;
    double clay_upper = MY_LARGE;
    for(int i=0;i<clay[1].size();i++){
        double z = msa->AtomByGlobalIndex(clay[1][i]).xyz[2];
        if( z < clay_upper ){
            clay_upper = z;
            clay_upper_index = clay[1][i];
        }
    }
    // calculate the total molecular weight of all atoms between two plates
    double totalMW = 0.0;
    // create a lookup table for quick access of atomic weights (also records what atoms are in the system)
    map<int,double> globleindex_atomicweights;
    for(int iMol=4;iMol<ms->MoleculesCount();iMol++){ // iMol0,1,2,3 are lower, upper plates and counter ions
        for(int iAtom=0;iAtom<(*ms)[iMol].AtomsCount();iAtom++){
            string ele = (*ms)[iMol][iAtom].element;
            double mw = PeriodicTable::AtomicWeight(ele);
            totalMW += mw;
            int globalindex = this->msa->GlobalIndexOfAtom((*ms)[iMol][iAtom]);
            globleindex_atomicweights[globalindex] = mw;
        }
    }
    // Calculate the average spacing.
    // Old way: simply minus a fixed space.
    /*
    double total_spacing = 0.0;
    for(int iFrame=0;iFrame<traj->NFrames();iFrame++){
        clay_upper = (*traj)[iFrame].X()[clay_upper_index][2];
        double spacing = clay_upper - clay_lower;
        total_spacing+=spacing;
    }
    total_spacing /= (traj->NFrames());
    total_spacing -= 6.0;  // subtract for untouchable spaces??
    // Calculate average density.
    double vol = ms->boundary.GetU()[0] * ms->boundary.GetV()[1] * total_spacing; // in A^3
    double rho = totalMW / vol / 6.022E23 * 1E24; // density in g/cm^3
    */
    // New way: Strict method
    double rho = 0.0;
    for(int iFrame=0;iFrame<traj->NFrames();iFrame++){
        clay_upper = (*traj)[iFrame].X()[clay_upper_index][2];
        double touchable_dist = 3.0;
        double lower = clay_lower+touchable_dist;
        double upper = clay_upper-touchable_dist;
        double spacing = upper-lower;
        double vol = ms->boundary.GetU()[0] * ms->boundary.GetV()[1] * spacing; // in A^3
        // atoms in between [lower, upper] are considered
        // calculate the total molecular weight of all atoms between two plates
        double totalMW = 0.0;
        for(auto &atomIndex_mw : globleindex_atomicweights){
            double z = (*traj)[iFrame].X()[atomIndex_mw.first][2];
            if(z<lower or z>upper)
                continue;
            totalMW += atomIndex_mw.second;
        }
        double this_rho = totalMW / vol / 6.022E23 * 1E24; // density in g/cm^3
        rho += this_rho;
    }
    rho /= (traj->NFrames());

    // output
    ofstream ofs;
    ofs.open(get_output_dir() / "RHO.csv");
    ofs << "# rho in g/cm^3 "<<endl;
    ofs<<rho<<endl;
}

/* Chain lengths distribution of all chains in all frames
 * the auxiliary function __chain_length__ calculates the length of chain at a certain frame */
double __chain_length__(TrajectoryFrame& frame,vector<int> atoms,double U,double V,double W){
    int nAtoms = atoms.size();
    double maxDist = 0.0;
    for(int i=0;i<nAtoms;i++){
        for(int j=i+1;j<nAtoms;j++){
            XYZ a1 = frame.X()[atoms[i]];
            XYZ a2 = frame.X()[atoms[j]];
            double d = __dist_in_PBC__(a1,a2,U,V,W);
            maxDist = max(maxDist,d);
        }
    }
    return maxDist;
}
void Analyzer::ChainLengths(double range_in_nm, double binSize_in_nm){
    vector<vector<int>> polymers;
    this->locate_polymers(polymers);
    int nChains = polymers.size(); // number of chains, should be 9
    int nBins = (int)ceil(range_in_nm/binSize_in_nm);
    vector<double> bins(nBins,0.0); // the bins to record length distributions
    double U = ms->boundary.GetU()[0];
    double V = ms->boundary.GetV()[1];
    double W = ms->boundary.GetW()[2];
    for(int iFrame=0;iFrame<traj->NFrames();iFrame++){
        TrajectoryFrame &f = (*traj)[iFrame];
        for(int iChain=0;iChain<nChains;iChain++){
            vector<int> &atoms = polymers[iChain];
            double d = __chain_length__(f,atoms,U,V,W) * 0.1; // A to nm.
            int iBin = (int)floor(d/binSize_in_nm);
            bins[iBin] += 1.0;
        }
    }
    for(int i=0;i<nBins;i++){
        bins[i] /= (traj->NFrames() * nChains);
    }
    // output
    ofstream ofs;
    ofs.open(get_output_dir() / "ChainLengths.csv");
    ofs << "#, NRows (dist), NColumns"<<endl;
    ofs<<nBins<<","<<2<<endl;
    for (int i = 0; i < nBins; i++)
        ofs<<(0.5+i)*binSize_in_nm<<", "<<bins[i]<<", "<<endl;
}

void __TestLocateModule__(string path,int max_workers){
    Analyzer ana(path,max_workers);
    vector<vector<int>> indexes;
    indexes.clear();ana.Locate(indexes,COO);
    ana.Locate_Demonstrate(indexes,"test_coo.mol2");
    indexes.clear();ana.Locate(indexes,AMMONIUM);
    ana.Locate_Demonstrate(indexes,"test_N+.mol2");
    indexes.clear();ana.Locate(indexes,CONH2);
    ana.Locate_Demonstrate(indexes,"test_conh2.mol2");
}

void MainAsAnalyze(string path,int max_workers){
    cout<<"Working on folder "<<path<<endl;
    Analyzer ana(path,max_workers);
    cout<<"PolymerZ..."<<flush;
    ana.PolymerZPositions();
    cout<<"ZDistrib..."<<flush;
    ana.ParticleZDistribution();
    cout<<"RDF..."<<flush;
    ana.RDF();
    cout<<"MSD..."<<flush;
    ana.MSD();
    cout<<"Rho..."<<flush;
    ana.Density();
    cout<<"ChainLengths..."<<flush;
    ana.ChainLengths();
    cout<<"All done."<<endl;
}

int main(int argc,char* argv[]){
    string path;
    int max_workers;
    if(argc<2)
        path = "./";
    else
        path = string(argv[1]);
    max_workers = argc<3 ? -1 : atoi(argv[2]);

    MainAsAnalyze(path,max_workers);
}
