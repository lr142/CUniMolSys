#ifndef CUNIMOLSYS_TRAJECTORY_H
#define CUNIMOLSYS_TRAJECTORY_H
#include "universalmolecularsystem.h"
#include <thread>
#include <mutex>
#include <future>
#include <random>
#include <set>
#include <deque>

/* The contents and structure of a specific LAMMPS trajectory file */
struct TrajFile{
    TrajFile() = default;
    ~TrajFile() = default;
    TrajFile(const TrajFile &other) = delete; // Prevent accidentally copied this huge object.
    TrajFile& operator=(const TrajFile &other) = delete;
    std::vector<string> lines; // All lines in the line. This is a huge data structure.
    string filename;
    vector<int> timesteps; // timestep of each frame; size of thie vector is # of frames in file.
    vector<int> startlines; // start line no of each frame, for quick access.
};
struct KeywordsColumnPos{
    /* The column of various fields in the LAMMPS trajectory file. Note that each frame may have
     * different fields and positions (if the dump file is writen by multiple different dump commands
     * in the LAMMPS script file).
     * The names of these variables are set to match those in LAMMPS dump commands:
     * possible attributes = id, mol, proc, procp1, type, element, mass,
        x, y, z, xs, ys, zs, xu, yu, zu,
        xsu, ysu, zsu, ix, iy, iz,
        vx, vy, vz, fx, fy, fz,
        q, mux, muy, muz, mu,
        radius, diameter, omegax, omegay, omegaz,
        angmomx, angmomy, angmomz, tqx, tqy, tqz,
        c_ID, c_ID[N], f_ID, f_ID[N], v_name
     * If a keyword is not present, its corresponding value is set to -1. Otherwise it's set to the column
     * which the keyword appears in each line of an atom.
     * We have could set a map<string,int> to store these number more elegantly, but they are written as
     *  variables for efficiency (saving the time to look up the map for millions of times )*/
    int id,mol,type,vx,vy,vz,fx,fy,fz,ix,iy,iz;
    int x,y,z,xu,yu,zu,xs,ys,zs,xsu,ysu,zsu;
    /* read a line in LAMMPS trajectory file like this:
     * ITEM: ATOMS id mol type x y z vx vy vz
     * and find the column number of each keyword. In the example above, the column of "id" is 0
     * and column of x is 3. */
    void FindColumnPos(string line);
};

/* Trajectory means a MolecularSystem with associated trajectories.
 * We use composite instead of inheritance since there may not be a situation where the user needs a MolecularSystem
 * while we only have a Trajectory at hand.
 * This class (in some sense) implements the Decorator design pattern */
class Trajectory{
public:
    Trajectory(MolecularSystem &ms);
    ~Trajectory();
    void _thread_main_(int iThread,int iFrameStart, int iFrameEnd,std::promise<bool> ret);
    void _testMultiThread();

    /* Read a trajectory file 'filename'. The caller can call Read() multiple times to read multiple files. The trajectory
     * will be appended to the last.
     * max_workers is the number of CPUs to use. <=0 means let the program itself to determine based on available resources
     * maxFrames is the maximum number of frames to actually read. default is unlimited
     * if removeDup==true, when a readly new frame has the same timestep as a previous frame (read from last
     *   call of this function), the old frame will be discarded. This is useful when a MD simulation writes multiple
     *   trajectory files, and sometimes the 1st frame of the 2nd file is exactly the last frame of the 1st file.
     * if certainFrames is non-empty, only the frames specified in certainFrames (by thier timestep values)  will be
     *   read. If the frame ts specified in [certainFrames] are not present in the actual file, this frame will be simply
     *   ignored. This is useful when the caller does not know how many frames are in the file.
     * This function returns the number of frames actually.
     * */
    int Read(string filename, int max_workers=-1, int maxFrames=(int)MY_LARGE,
             bool removeDup=true, std::set<int> certainFrames=set<int>());
    bool DiscardFrame(int iFrame);
protected:
    MolecularSystem &ms_;
    /* number of atoms in the system. Assume all frames have the same number of atoms. This means GCMC is not directly
    supported. However, enough space can be reserved in the beginning for atoms coming into the system at a later time. */
    int nAtoms_;
    int nFrames_; // number of frames
    double timestep_in_fs; // timestep in fs (Same as the convention of 'units real' in LAMMPS)

    vector<int> ts_; // size() == nFrames_; recorded timesteps of each frame.
    /* The following data structures stores the core info. Each of which is a vector of length nFrames. e.g. x_[i] records
     * the coordinates of all atoms in frame i. x_[i] itself is an array of struct XYZ with length = nAtoms.
     * Note: velocities/forces/image_flags may not be available in all trajectories.
     * But these vectors are **ALWAYS** having lengths of nFrames.
     * If in some frames, the corresponding info is absent, say i_ is not present in frame 10, then
     * i_[9] == nullptr;
     * Currently using double to store these data. If memory is an issue, recompile with using XYZ =  XYZ_T_<float> in utility.h;
     * */
    vector<XYZ*> x_; // coordinates
    vector<XYZ*> v_; // velocities
    vector<XYZ*> f_; // forces;
    vector<XYZ_T_<int>*> i_; // periodic image flags (ix,iy,iz in LAMMPS).
    /* x_,v_,f_,i_ are highly isomorphic. I wrote this macro to carry out same operations on all vectors instead of writing
     * these four vector again and again. If in the future another vector is added, e.g. atomic charges (which is changing)
     * add [vector<double*> charges_;] above, and add [opers(charges_,(arg))] in the macro below;
     * Similarly, ts_ records info about the system as every trajectory frame. Another macro is defined to perform
     * operations on this kind of vector for possibility of future extension (e.g. system size in PVT sims);
     * */

#define OPERATION_FOR_ALL_PERATOM_VECTORS(opers,arg) opers(x_,(arg)); opers(v_,(arg)); opers(f_,(arg)); opers(i_,(arg));
#define OPERATION_FOR_ALL_SYSTEM_VECTORS(opers,arg) opers(ts_,(arg));
    /* For your convenience, below is a list of places you need to modify when new properties are to be included:
     * 1. Declare this property in Trajectory::vector<T> newProp;
     * 2. If this is a system-wise property (one value for each frame), add its name in OPERATION_FOR_ALL_SYSTEM_VECTORS macro
     * 3. If this is a per-atom property (like XYZ):
     * 3.1 Add its name in OPERATION_FOR_ALL_PERATOM_VECTORS macro and
     * 3.2 Declare its name in KeywordsColumnPos. And Modify KeywordsColumnPos::FindColumnPos(string);
     * 3.3 Modify Trajectory::createMemoryForFrame(), when should the memory be created?
     * 4. Modify Trajectory::read_one_frame(), how should the new property be read and recorded?
     * */

protected:
    void createMemoryForFrame(int iFrame,KeywordsColumnPos &kcp);
    void destroyMemoryForFrame(int iFrame);
    void destroyAll();

    /* read file contents into trajFile */
    void read_file_step0(string filename);
    /* Given that all lines of a lammpstrj file is read into this->trajFile, this function
     * skims the contents of the file and find out how many frames are there and
     * what are the timesteps (trajFile.timesteps) and starting line numbers (trajFile.startlines) of each frame.
     * This function */
    void read_preparation_step1_find_frames_in_file();
    /* Based on the given maxFrames and/or certainFrames parameters, determine which frames will be actually read.
     * This function will modify (trajFile.timesteps) and (trajFile.startlines) to reflect the result */
    void read_preparation_step2_find_actually_read_frames(int maxFrames, set<int> &certainFrames);
    /* This function will remove duplicate frames in the existing data. By duplication we mean two frames have the
     * same timestep value. The new frames are given in (trajFile.timesteps). If duplication occurs, the old data
     * will be deleted */
    void read_preparation_step3_remove_duplication(bool removeDup);

    void read_one_frame(int iFrame, int iFrameInTrajFile);
    TrajFile trajFile;
    std::default_random_engine e;
};



#endif //CUNIMOLSYS_TRAJECTORY_H

