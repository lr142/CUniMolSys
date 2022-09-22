#ifndef CUNIMOLSYS_TRAJECTORY_H
#define CUNIMOLSYS_TRAJECTORY_H
#include "universalmolecularsystem.h"
#include <thread>
#include <mutex>
#include <future>
#include <random>
#include <set>
#include <deque>
#include <memory>
using std::vector;
using std::shared_ptr;

/* The contents and structure of a specific LAMMPS trajectory file */
class TrajFile{
public:
    TrajFile() = default;
    ~TrajFile() = default;
    TrajFile(const TrajFile &other) = delete; // Prevent accidentally copied this huge object.
    TrajFile& operator=(const TrajFile &other) = delete;
    int NFrames();
    void AddFrame(int nAtoms,int ts);
    inline vector<string>& operator[] (int iFrame) { return *(lines[iFrame]);}
    inline int NAtoms(int iFrame){return nAtoms[iFrame];}
    inline int TimeStep(int iFrame){ return timesteps[iFrame];}
    void Clear();
    string filename;
protected:
    /* All lines in the line.lines for each frame are stored in a vector<string>, and different frames are stored
     * as a vector of shared_ptr of vector<string>. I tried to store all lines in a huge vector<string> but that
     * was a bad idea */
    vector<shared_ptr<vector<string>>> lines;
    vector<int> timesteps; // timestep of each frame; size of this vector is # of frames in file.
    vector<int> nAtoms; // number of atoms of each frame in this file.
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

/* A frame in the trajectory */
class TrajectoryFrame{
public:
    TrajectoryFrame(int nAtoms);
    ~TrajectoryFrame();
    TrajectoryFrame(const TrajectoryFrame& tra) = delete; // prevent accidental copy
    TrajectoryFrame& operator=(const TrajectoryFrame& tra) = delete;
    /* Read a frame from the file. the index of frame in the file is iFrameInTrajFile. if createMemory is true, the
     * function will create memory itself based on the read atoms. If not, the memory is already ready. */
    void Read(TrajFile &trajFile, int iFrameInTrajFile);
    /* number of atoms in this frame. This number may be smaller (if user dump a specific group instead of all)
 * or larger (when running GCMC) then number of atoms in the system. */
    int nAtoms_;
    int ts_; // timestep of this frame.
    /* */
    int* s_; // global serials of every atom in this frame. ('id' reported in dump file)
    XYZ* x_; // coordinates
    XYZ* v_; // velocities
    XYZ* f_; // forces;
    XYZ_T_<int>* i_; // periodic image flags (ix,iy,iz in LAMMPS).
    friend class Trajectory;
protected:
    /* Create memory, allocate space as decribed in kcp. Before called this, nAtoms_ must be correctly initialized */
    void createMemory(KeywordsColumnPos &kcp);
    void destroyMemory();
    void sort_atoms();




};

/* Trajectory means a MolecularSystem with associated trajectories.
 * We use composite instead of inheritance since there may not be a situation where the user needs a MolecularSystem
 * while we only have a Trajectory at hand.
 * This class (in some sense) implements the Decorator design pattern */
class Trajectory{
public:
    Trajectory(MolecularSystem &ms);
    ~Trajectory();
    void _thread_main_(int iThread,int iFrameStart, int iFrameEnd, int NAtoms);
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
     * This function returns the number of actually read frames.
     * */
    int Read(string filename, int max_workers=-1, int maxFrames=(int)MY_LARGE,
             bool removeDup=true, std::set<int> certainFrames=set<int>());
    /* Discard a specific frame */
    bool DiscardFrame(int iFrame);
    /* Discard all frames, release the memory, and get ready to another trajctory file */
    void Clear();
    int NFrames();

    /* Returns a copy of the molsys at iFrame with all coordinates updated to this frame.
     * Atoms not involved in the trajectory will keep their original coords. If only_in_traj
     * is true, only those atoms reported in the trajectory will be reported */
    shared_ptr<MolecularSystem> UpdateCoordsAtFrame(int iFrame, bool only_in_traj=false);
    void ShowTrajectory(string filename, bool showOriginal=true);
    inline TrajectoryFrame & operator[](int i) {return *(frames_[i]);}
protected:
    MolecularSystem &ms_;

    /* all frames. Tried to implement this as vector<TrajectoryFrame>. But as I tried to resize it, the previous frames
     * are haveing memeory issues (previously allocated memory lost). */
    vector<shared_ptr<TrajectoryFrame>> frames_;

    /* Read all lines of a lammpstrj file into this->trajFile.
     * skims the contents of the file and find out how many frames are there and
     * what are the timesteps (trajFile.timesteps) and starting line numbers (trajFile.startlines) of each frame.
     * This function */
    void read_preparation_step1_find_frames_in_file(string filename);
    /* Based on the given maxFrames and/or certainFrames parameters, determine which frames will be actually read.
     * This function will modify (trajFile.timesteps) and (trajFile.startlines) to reflect the result */
    void read_preparation_step2_find_actually_read_frames(int maxFrames, set<int> &certainFrames);
    /* This function will remove duplicate frames in the existing data. By 'duplication' we mean two frames have the
     * same timestep value. The new frames are given in (trajFile.timesteps). If duplication occurs, the old data
     * will be deleted */
    void read_preparation_step3_remove_duplication(bool removeDup);

    void read_initialize_step4_initialize_frame_vector(int &oldNFrames);
    void read_final_step5_in_parallel(int oldNFrames, int max_workers);

    TrajFile trajFile;
    std::default_random_engine e;

    void __test_thread_main__(int iThread,int &iFrame,int nFrames, int nAtoms);

    void __read_frame_thread_main__(int iThread, int NThreads, int oldNFrames);

    std::mutex mux;
};



#endif //CUNIMOLSYS_TRAJECTORY_H

