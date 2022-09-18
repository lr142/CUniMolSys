#ifndef CUNIMOLSYS_MOLECULARMANIPULATOR_H
#define CUNIMOLSYS_MOLECULARMANIPULATOR_H
#include "utility.h"
#include "universalmolecularsystem.h"

/* This function reorganizes atoms/molecules in a MolecularSystem ms
 * scheme is a vector with length = atoms count in ms, i-th atom in the old molecular system
 * will belong to scheme[i]-th molecule in the new molecular system.
 * The user should make sure 0<=scheme[i]<molecules count in the new system.
 * There may be cases where some molecules in the new system contains no atom. eg. scheme = [0,2]
 * The original MolecularSystem ms will be destroyed, the result will be writing into ms.
 * If the caller wants to retain the original copy, they should call ms.Deepcopy() to save a copy
 * before calling this function.
 */
void MolSysReorganize(MolecularSystem &ms, vector<int> &scheme);

/* Randomly split the molecular system into N molecules. Used mainly for testing. */
void MolSysRandomSplit(MolecularSystem &ms, int N);

/* Combine multiple molecules in a molecular system into a single molecule */
void MolSysReduceToSingleMolecule(MolecularSystem &ms);
#endif //CUNIMOLSYS_MOLECULARMANIPULATOR_H
