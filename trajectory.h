#ifndef CUNIMOLSYS_TRAJECTORY_H
#define CUNIMOLSYS_TRAJECTORY_H
#include "universalmolecularsystem.h"

/* Trajectory means a MolecularSystem with associated trajectories.
 * We use composite instead of inheritance since there may not be a situation where the user needs a MolecularSystem
 * while we only have a Trajectory at hand.
 * This class (in some sense) implements the Decorator design pattern */

class Trajectory{
public:
    Trajectory(MolecularSystem &ms);
    ~Trajectory();
protected:
    MolecularSystem &ms_;
    /* number of atoms in the system. Assume all frames have the same number of atoms. This means GCMC is not directly
    supported. However, enough space can be reserved in the beginning for atoms coming into the system at a later time. */
    int nAtoms_;
    int nFrames_; // number of frames
    vector<int> timesteps_each_frame_; // size() == nFrames_;
    double timestep_in_fs; // timestep in fs (Same as the convention of 'units real' in LAMMPS)

    /* The following data structures stores the core info. Each of which is a vector of length nFrames. e.g. x_[i] records
     * the coordinates of all atoms in frame i. x_[i] itself is an array of struct XYZ with length = nAtoms.
     * Note: velocities/forces/image_flags may not be available in all trajectories. In those cases, they are simply
     * set to empty vectors.
     * Current using double to store these data. If memory is an issue, recompile with XYZ_T_<float>;
     * */
    vector<XYZ*> x_; // coordinates
    vector<XYZ*> v_; // velocities
    vector<XYZ*> f_; // forces;
    vector<XYZ_T_<int>*> image_flags_; // periodic image flags (ix,iy,iz in LAMMPS).
protected:
    void createMemoryForFrame(int iFrame);
    void destroyMemoryForFrame(int iFrame);
    void destroyAll();
};
#endif //CUNIMOLSYS_TRAJECTORY_H

