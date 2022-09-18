#include "molecularmanipulator.h"
#include <random>
#include <stack>
using namespace std;
void MolSysReorganize(MolecularSystem &ms, vector<int>& scheme){
    MolecularSystemAccessor oldMolAccessor(ms);
    // validity check
    int nAtoms = ms.AtomsCount();
    if(nAtoms == 0) // Empty System
        return;
    if(scheme.size()!=nAtoms)
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
    for(int iMol=0;iMol<newMolsContainedAtomIndexes.size();iMol++){
        for(int iAtom=0;iAtom<newMolsContainedAtomIndexes[iMol].size();iAtom++){
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
