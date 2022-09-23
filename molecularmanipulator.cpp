#include "molecularmanipulator.h"
#include "opener.h"
#include "bonddetector.h"
#include "utility.h"
#include <random>
#include <stack>
#include <string>
#include <sstream>
using namespace std;
void MolSysReorganize(MolecularSystem &ms, vector<int>& scheme){
    MolecularSystemAccessor oldMolAccessor(ms);
    // validity check
    int nAtoms = ms.AtomsCount();
    if(nAtoms == 0) // Empty System
        return;
    if((int)scheme.size()!=nAtoms)
        ERROR("Length of scheme must be equal to atoms count in the MolecularSystem");
    int newMolsCount = 0;
    for(auto item:scheme){
        if(item<0)
            ERROR("scheme[i] must not be negative");
        newMolsCount = max(newMolsCount,item+1);
    }
    // Create a temp new MolecularSystem newMS
    MolecularSystem newMS;
    newMS.molecules.resize(newMolsCount);
    for(int iMol=0; iMol < newMolsCount; iMol++){
        newMS.molecules[iMol] = make_shared<Molecule>();
    }
    //  Copy atoms.
    vector<vector<int>> newMolsContainedAtomIndexes(newMolsCount);
    for(int iAtom=0;iAtom<nAtoms;iAtom++){
        int iBelongsToWhichNewMol = scheme[iAtom];
        auto locationInOldMol = oldMolAccessor.MolAndLocalIndexOfAtom(iAtom);
        int oldMolIndex = locationInOldMol.first;
        int oldLocalIndex = locationInOldMol.second;
        newMS[iBelongsToWhichNewMol].atoms.push_back(ms[oldMolIndex].atoms[oldLocalIndex]);
        // Meanwhile, we also need to keep a record of the atoms each new molecule contains
        newMolsContainedAtomIndexes[iBelongsToWhichNewMol].push_back(iAtom);
    }
    // Next we build two important maps, globalIndexOldToNew and globalIndexNewToOld
    vector<int> globalIndexOldToNew(scheme.size());
    vector<int> globalIndexNewToOld(scheme.size());
    int newIndexCounter = 0;
    for(unsigned int iMol=0;iMol<newMolsContainedAtomIndexes.size();iMol++){
        for(unsigned iAtom=0;iAtom<newMolsContainedAtomIndexes[iMol].size();iAtom++){
            int oldIndex = newMolsContainedAtomIndexes[iMol][iAtom];
            globalIndexOldToNew[oldIndex] = newIndexCounter;
            globalIndexNewToOld[newIndexCounter] = oldIndex;
            ++newIndexCounter;
        }
    }
    // Renumber the newMS and create its Accessor
    newMS.RenumberAtomSerials();
    MolecularSystemAccessor newMolAccessor(newMS);
    // Next we take care of the bonds. Read all bonds from oldMolAccessor
    auto & bondedMap = oldMolAccessor.GetBondedMap();
    for(int iOldFromAtom=0;iOldFromAtom<nAtoms;iOldFromAtom++){
        for(auto& item:bondedMap[iOldFromAtom]){
            int iOldToAtom = item.first;
            // store only iOldFromAtom < iOldToAtom
            if(iOldFromAtom >= iOldToAtom)
                continue;
            shared_ptr<Bond> pbond = item.second;
            // For the following variables, [0] and [1] refers to from and to.
            int iOldAtoms[2] = {iOldFromAtom,iOldToAtom}; // globalIndexes of atoms in old system
            int iNewMols[2]; // mol indexes of two involved atoms;
            string localSerials[2]; // local indexes of two involved atoms;
            string globalSerials[2]; // globalSerial of two involved atoms;
            for(int i=0;i<2;i++){
                int iNewAtoms = globalIndexOldToNew[iOldAtoms[i]];
                auto iPairs = newMolAccessor.MolAndLocalIndexOfAtom(iNewAtoms);
                iNewMols[i] = iPairs.first;
                int iLocalIndex = iPairs.second;
                localSerials[i] = newMS[iNewMols[i]][iLocalIndex].serial;
                globalSerials[i] = newMS[iNewMols[i]][iLocalIndex].globalSerial;
            }
            Bond b;
            b.length = pbond->length;
            b.type = pbond->type;
            // in-molecules bonds;
            if(iNewMols[0]==iNewMols[1]){
                b.atom1 = localSerials[0];
                b.atom2 = localSerials[1];
                newMS[iNewMols[0]].bonds.push_back(make_shared<Bond>(b));
            }else{
                b.atom1 = globalSerials[0];
                b.atom2 = globalSerials[1];
                newMS.interMolecularBonds.push_back(make_shared<Bond>(b));
            }
        }
    }
    // Simply replace the molecules array and inter-molecular bonds array with the new ones
    // will be enough. Thanks to shared_ptrs.
    // BTW, the trajectory(if there is) will become invalid since atom orders are different.
    // but this is not handed by this function.
    ms.molecules = newMS.molecules;
    ms.interMolecularBonds = newMS.interMolecularBonds;
    ms.RenumberAtomSerials();
}

void MolSysRandomSplit(MolecularSystem &ms, int N) {
    int nAtoms = ms.AtomsCount();
    if(nAtoms==0)
        return;
    if(N<1)
        ERROR("N must >= 1");
    default_random_engine e;
    e.seed(time(0));
    uniform_int_distribution<int> uid(0,N-1);
    vector<int> scheme(nAtoms);
    for_each(scheme.begin(),scheme.end(),
             [&uid,&e](decltype(*scheme.begin()) &it){it=uid(e);});
    // This function will guarantee that the maximum molecular index == N-1 by manually setting the last element
    scheme[scheme.size()-1] = N-1;
    MolSysReorganize(ms,scheme);
}

void MolSysReduceToSingleMolecule(MolecularSystem &ms) {
    vector<int> scheme(ms.AtomsCount(),0);
    MolSysReorganize(ms,scheme);
}

void MolSysSplitByConnectivity(MolecularSystem &ms) {
    // This is a standard graph exploration
    // We will use Depth-First-Search (DFS) for simplicity
    int nAtoms = ms.AtomsCount();
    if(nAtoms==0)
        return;
    MolecularSystemAccessor msa(ms);

    vector<bool> visited(nAtoms,false);
    vector<int> parent(nAtoms); // parent[i] == i means this is a root node (or unvisited)
    for(int i=0;i<nAtoms;i++)
        parent[i] = i; // parent == itself means a root node
    for(int root=0;root<nAtoms;root++){
        if(visited[root])
            continue;
        visited[root] = true;
        stack<int> to_visit;
        to_visit.push(root);
        while(!to_visit.empty()){
            int node = to_visit.top();
            to_visit.pop();
            for(auto &item:msa.GetBondedMap()[node]){
                int iToAtom = item.first;
                if(visited[iToAtom])
                    continue;
                else{
                    visited[iToAtom] = true;
                    parent[iToAtom] = node;
                    to_visit.push(iToAtom);
                }
            }
        }

    }
    // Find the root of every atom
    for(int i=0;i<nAtoms;i++){
        while(parent[parent[i]] != parent[i])
            parent[i] = parent[parent[i]];
    }
    // Find how many roots (trees) are there
    int treeCount = 0;
    map<int,int> rootToTreeMap;
    for(int i=0;i<nAtoms;i++){
        int root= parent[i];
        if(rootToTreeMap.find(root) == rootToTreeMap.end())
            rootToTreeMap.insert(make_pair(root,treeCount++));
    }
    // Finally generate the corresponding vector
    vector<int> scheme(nAtoms);
    for(int i=0;i<nAtoms;i++){
        scheme[i] = rootToTreeMap[parent[i]];
    }
    MolSysReorganize(ms,scheme);
}

void MolSysExtend(MolecularSystem &dest, MolecularSystem &src){
    int nAtomsSysInDest = dest.AtomsCount();
    int nInterBondsInDest = dest.interMolecularBonds.size();

    // Just copy the shared_ptrs, not really the Atoms, bonds. etc.
    dest.molecules.insert(dest.molecules.end(),src.molecules.begin(),src.molecules.end());
    // Intermolecular bonds of src mol should be treated differently
    MolecularSystemAccessor msaOld(src);
    dest.interMolecularBonds.insert(dest.interMolecularBonds.end(),
                                    src.interMolecularBonds.begin(),
                                    src.interMolecularBonds.end());
    dest.RenumberAtomSerials();
    MolecularSystemAccessor msaNew(dest);
    for(unsigned int iBond=nInterBondsInDest;iBond<dest.interMolecularBonds.size();iBond++){
        auto pNewBond = dest.interMolecularBonds[iBond];
        string oldGlobalSerials[2];
        string newGlobalSerials[2];
        for(unsigned int i=0;i<2;i++){
            int oldGlobalIndex = msaOld.GlobalIndexOfAtom(oldGlobalSerials[i]);
            int newGlobalIndex = oldGlobalIndex + nAtomsSysInDest;
            newGlobalSerials[i] = msaNew.AtomByGlobalIndex(newGlobalIndex).globalSerial;
        }
        pNewBond->atom1 = newGlobalSerials[0];
        pNewBond->atom2 = newGlobalSerials[1];
    }
}

void MolSysDuplicatePeriodically(MolecularSystem &ms, int ix, int iy, int iz, bool set_bound){
    if(! ms.Periodic())
        ERROR("System not periodic.");
    if(ix<1 or iy<1 or iz<1)
        ERROR("Image count must be >=1 ");
    MolecularSystem original = ms.DeepCopy();

    int nTotal = ix*iy*iz;
    int iCount = 0;
    bool requiresProgressBar = false;
    if(nTotal*ms.AtomsCount() > MY_LARGE_SYSTEM)
        requiresProgressBar = true;  // Requires Progress Bar if over 0.5 million atoms;


    for(int i=0;i<ix;i++){
        for(int j=0;j<iy;j++){
            for(int k=0;k<iz;k++){
                if(i==0 and j==0 and k==0)
                    continue;
                MolecularSystem copy = original.DeepCopy();
                XYZ offset = original.boundary.GetU()*i +
                        original.boundary.GetV()*j +
                        original.boundary.GetW()*k;
                copy.Translate(offset);
                MolSysExtend(ms,copy);
                if(requiresProgressBar)
                    ProgressBar(1.0*iCount/nTotal);
                ++iCount;
            }
        }
    }
    if(set_bound){
        auto u = ms.boundary.GetU()*ix;
        auto v = ms.boundary.GetV()*iy;
        auto w = ms.boundary.GetW()*iz;
        ms.boundary.SetUVW(u,v,w);
    }
    ms.RenumberAtomSerials();
    if(requiresProgressBar) {
        ProgressBar(1.0);
    }
}

void MolSysSubsystemByMask(MolecularSystem &ms,vector<bool> mask){
    if((int)mask.size() != ms.AtomsCount())
        ERROR("vector mask must have same length as the number of atoms in MolecularSystem");
    vector<int> scheme(mask.size());
    for(unsigned int i=0;i<scheme.size();i++){
        scheme[i] = mask[i]? 0 : 1;
    }
    MolSysReorganize(ms,scheme);
    // Keep only the first molecule.
    ms.molecules.erase(ms.molecules.begin()+1,ms.molecules.end());
    // Discard all inter-molecular bonds
    ms.interMolecularBonds.clear();
    ms.RenumberAtomSerials();
}

void MolSysSubSystemBySerials(MolecularSystem &ms,set<string> &globalSerials){
    vector<bool> mask(ms.AtomsCount());
    int indexCounter = 0;
    for(int i=0;i<ms.MoleculesCount();i++){
        for(int j=0;j<ms[i].AtomsCount();j++){
            if(globalSerials.contains(ms[i][j].globalSerial))
                mask[indexCounter] = true;
            else
                mask[indexCounter] = false;
            ++indexCounter;
        }
    }
    MolSysSubsystemByMask(ms,mask);
}

bool _add_this_water_(MolecularSystemAccessor &originalSys, NeighborList &nlist, Boundary &bound, Molecule &mol, int indexOfOxygen, double minDistSquared){
    XYZ oxygen_pos = mol[indexOfOxygen].xyz;
    XYZ_DTYPE lohi[3][2];
    bound.GetLoHi(lohi);
    // Check if the water is out of the bound
    for(int i=0;i<3;i++){
        if(oxygen_pos[i] < lohi[i][0] or oxygen_pos[i]>lohi[i][1])
            return false;
    }
    // Clash Detection
    auto neighbors = nlist.GetNeighborsFromCachedLists(oxygen_pos);
    if(neighbors== nullptr)
        return true;

    for(auto iter=neighbors->begin();iter!=neighbors->end();++iter){
        Atom &candidateAtom = originalSys.AtomByGlobalIndex(iter->index);
        /* Must call nlist.DistanceSquared() instead of comparing these two coords directly due to
         * 1. In PBC, the waters coordinates must be wrapped in the cell to detect clash with system atoms in another image
         * 2. In Non-PBC, the coordinates from nlist->iter->xyz are shifted away from the original by skin_size/2 */
        double distSqrt = nlist.DistanceSquared(oxygen_pos,candidateAtom.xyz);
        if(distSqrt < minDistSquared)
            return false;
    }
    return true;
}

bool MolSysSolvate(MolecularSystem &ms, Boundary bound, WaterType waterType, double minDistance){
    if(!bound.Orthogonal())
        ERROR("Supports only orthogonal region");
    string path = DATAFILESPATH+"/Structures/";
    if(waterType==WaterType::SPC)
        path += "spc_water_5nm_chunk.mol2";
    else if(waterType==WaterType::TIP4P)
        path += "tip4p_water_5nm_chunk.mol2";
    else
        ERROR("Unknown water type"+ to_string(waterType));

    MolecularSystem waterChunk;
    QuickOpen(waterChunk,path);
    waterChunk.ClearBonds();
    if(waterChunk.AtomsCount()<3)
        ERROR("No waters read from the water chunk file");
    // Erase atom type info in the .mol2 file, auto recognize later
    for(int i=0;i<waterChunk.MoleculesCount();i++){
        for(int j=0;j<waterChunk[i].AtomsCount();j++){
            waterChunk[i][j].type = "";
        }
    }
    // Calc the fill size
    double chunkSize = waterChunk.boundary.GetU()[0]; // The chunk should be a cubic, usally of size 50 Angs.
    // Align the origin of the waterchunk to the desired origin. The origin of the water chunk is usually [0,0,0], but no
    // harm to just write this.
    waterChunk.Translate(bound.GetOrigin()-waterChunk.boundary.GetOrigin());

    XYZ fillSize = { bound.GetU()[0], bound.GetV()[1], bound.GetW()[2] };
    int periodicImages[3];
    for(int i=0;i<3;i++){
        periodicImages[i] = (int)ceil(fillSize[i]/chunkSize);
        if(periodicImages[i]<=0)
            ERROR("Fill region size can't be 0 or negative");
        else if(fillSize[i]>=1000) // Warning if region > 1000 Angstrom
            WARNING("Fill region too large : "+ to_string(fillSize[i])+" Angstroms, are you sure?");
        else{}
    }
    // Prepare the chunk for filling. Duplicate Images:
    ostringstream msg;
    msg<<"Creating water chunk of size = "<<bound.GetU()[0]/10.0<<"*"<<bound.GetV()[1]/10.0<<"*"<<bound.GetW()[2]/10.0<<" nm...";
    output(msg.str());
    MolSysDuplicatePeriodically(waterChunk,periodicImages[0],periodicImages[1],periodicImages[2],true);
    MolSysSplitByConnectivity(waterChunk); // Each water will occupy a molecule

    // Find the index of O atom in each molecule. All waters are the same so this is done only once.
    int indexOfOxygen;
    for(indexOfOxygen=0;indexOfOxygen<waterChunk[0].AtomsCount();indexOfOxygen++){
        if(StringToUpper(waterChunk[0][indexOfOxygen].element)=="O")
            break;
    }

    // Add water mols into the system
    NeighborList nlist(ms,minDistance*1.5); // gridSize = minDistance should be enough, *1.5 just for safety.
    nlist.BuildCachedNeighborList(); // Because we will use cached neighbor list, not the precise list

    MolecularSystem tempSys;
    MolecularSystemAccessor originalSysAccessor(ms);
    double minDistSquared = pow(minDistance,2);
    output("Adding Waters to System...");
    for(int i=0;i<waterChunk.MoleculesCount();i++){
        if(_add_this_water_(originalSysAccessor,nlist,bound,waterChunk[i],indexOfOxygen,minDistSquared)){
            tempSys.AddMolecule(waterChunk[i]);
        }
    }

    if(tempSys.MoleculesCount()==0) // No water added.
        return false;

    MolSysReduceToSingleMolecule(tempSys);
    ms.AddMolecule(tempSys[0]);
    ms.RenumberAtomSerials();
    output("Solvating the System Done.");
    return true;
}

