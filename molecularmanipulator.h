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

/* Reorganize a MolecularSystem into multiple molecules according to connectivity.
 * May increase or decrease the number of molecules */
void MolSysSplitByConnectivity(MolecularSystem &ms);

/* Extend the system by copied molecules, bonds, and inter-molecular bonds to dest
 * The contents of src molecules are used by dest (the shared_ptrs are copied)
 * Therefore if the user wants to continue using src after this function, the user
 * should pass a src.Deepcopy() instead of src into this function
 */
void MolSysExtend(MolecularSystem &dest, MolecularSystem &src);

/* Duplicate a system periodically.
 * ix,iy,iz must >= 1
 * Works only for periodic systems. Supports non-orthogonal system!
 * the boundary info of the ms will also be modified if set_bound flag is set to true */
void MolSysDuplicatePeriodically(MolecularSystem &ms, int ix, int iy, int iz, bool set_bound);

/* Make a subsystem by picking all atoms i which mask[i] == true.
 * Will put all picked-atoms in one molecule */
void MolSysSubsystemByMask(MolecularSystem &ms,vector<bool> mask);
/* Make a subsystem by picking serials */
void MolSysSubSystemBySerials(MolecularSystem &ms,set<string> &globalSerials);


/*
 *  Solvate a system by adding water molecules.
 *  bound is the region to add water. Supports Orthogonal region only
 *  watertype now supports  SPC and TIP4P. The coordinates of solvate molecules are pre-equilibrated water chunks
 *  of size 5nm*5nm*5n, located in DATAFILESPATH/Structures/
 *  minDistance is the minimal allowed distance between water O atom and any atom in the current system
 *  If successful, the MolecularySystem will has an additional molecule containing all waters. If failed, the ms will
 *  be untouched, and returns false.
 *  Note: this function will not check the periodicity of the system. The caller should make sure the ms is large enough
 *  to contain the region defined by bound.
  */
enum WaterType {TIP4P, SPC};
bool MolSysSolvate(MolecularSystem &ms, Boundary bound, WaterType waterType, double minDistance);

#endif //CUNIMOLSYS_MOLECULARMANIPULATOR_H

