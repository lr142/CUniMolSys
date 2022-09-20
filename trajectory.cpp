#include "trajectory.h"
#include "opener.h"
#include <future>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

// Auxiliary functions used in this module
void my_assert_true(bool statement){
    if(not (statement))
        throw exception();
}
bool double_equal(double a,double b){
    return fabs(a-b)<MY_SMALL;
}
#define MY_ASSERT_DOUBLE_EQUAL(a,b) if(not double_equal((a),(b))){cout<<iFrame<<" "<<iAtom<<" "<<(a)<<" "<<(b)<<endl;}
template <class T> void vector_free_space_at_pos(vector<T*> &vec, int iFrame, int nAtoms, bool missing_condition){
    if(missing_condition)
        vec[iFrame] = nullptr;
    else
        vec[iFrame] = new T[nAtoms];

}
template <class T> void vector_create_space_at_pos(vector<T*> &vec, int iFrame){
    if(vec[iFrame] != nullptr) {
        delete[] vec[iFrame];
        vec[iFrame] = nullptr;
    }
}
template <class T> void vector_erase_at(vector<T> &vec, int pos){
    if(vec.size()!=0) {
        my_assert_true(vec.size() > pos);
        vec.erase(vec.begin() + pos, vec.begin() + pos + 1);
    }
}
template <class T> void vector_resize(vector<T> &vec, int size){
    vec.resize(size);
}

void KeywordsColumnPos::FindColumnPos(string line) {
    // a line like this: "ITEM: ATOMS id mol type x y z vx vy vz"
    auto parts = StringSplit(line,' ');
    parts.erase(parts.begin(),parts.begin()+2);
    map<string,int> k2c; // key to column pos map
    for(int i=0;i<parts.size();i++){
        k2c[parts[i]] = i;
    }
#define FIND_POS(VAR) VAR = k2c.contains(#VAR) ? k2c[#VAR] : -1;
    FIND_POS(id) FIND_POS(mol) FIND_POS(type)
    FIND_POS(vx) FIND_POS(vy) FIND_POS(vz)
    FIND_POS(fx) FIND_POS(fy) FIND_POS(fz)
    FIND_POS(ix)  FIND_POS(iy)  FIND_POS(iz)
    FIND_POS(x)  FIND_POS(y)  FIND_POS(z)
    FIND_POS(xu)  FIND_POS(yu)  FIND_POS(zu)
    FIND_POS(xs)  FIND_POS(ys)  FIND_POS(zs)
    FIND_POS(xsu)  FIND_POS(ysu)  FIND_POS(zsu)
#undef FIND_POS
    // Special rules: if x not reported but xu (unwrapped), xs(scaled), xsu(scaled unwrapped)
    // Use that in the line of succession x <-- xu <-- xs <-- xsu
#define REPLACE_IF(VAR1, VAR2) if (VAR1==-1 and VAR2!=-1){VAR1=VAR2;}
    REPLACE_IF(xs,xsu)  REPLACE_IF(xu,xs)  REPLACE_IF(x,xu)
    REPLACE_IF(ys,ysu)  REPLACE_IF(yu,ys)  REPLACE_IF(y,yu)
    REPLACE_IF(zs,zsu)  REPLACE_IF(zu,zs)  REPLACE_IF(z,zu)
#undef REPLACE_IF
}

Trajectory::Trajectory(MolecularSystem &ms):ms_(ms),e(0){
    nAtoms_ = ms.AtomsCount();
    nFrames_ = 0;
}
Trajectory::~Trajectory(){
    destroyAll();
}

void Trajectory::createMemoryForFrame(int iFrame,KeywordsColumnPos &kcp) {
    bool cond;
    // Unfortunately, we can't write a macro to make this function simpler.
    cond = kcp.x<0 and kcp.y<0 and kcp.z<0;
    vector_free_space_at_pos(x_, iFrame, nAtoms_, cond);

    cond = kcp.vx<0 and kcp.vy<0 and kcp.vz<0;
    vector_free_space_at_pos(v_, iFrame, nAtoms_, cond);

    cond = kcp.fx<0 and kcp.fy<0 and kcp.fz<0;
    vector_free_space_at_pos(f_, iFrame, nAtoms_, cond);

    cond = kcp.ix<0 and kcp.iy<0 and kcp.iz<0;
    vector_free_space_at_pos(i_, iFrame, nAtoms_, cond);
}

void Trajectory::destroyMemoryForFrame(int iFrame){
    OPERATION_FOR_ALL_PERATOM_VECTORS(vector_create_space_at_pos, iFrame);
}
void Trajectory::destroyAll() {
    int nTotalFrame = this->nFrames_;
    for(int i=nTotalFrame-1;i>=0;i--)
        this->DiscardFrame(i);
    my_assert_true(this->nFrames_ == 0);
}

bool Trajectory::DiscardFrame(int iFrame){
    if(iFrame<0 or iFrame>=nFrames_)
        return false;
    destroyMemoryForFrame(iFrame);
    OPERATION_FOR_ALL_PERATOM_VECTORS(vector_erase_at, iFrame);
    OPERATION_FOR_ALL_SYSTEM_VECTORS(vector_erase_at, iFrame);
    --nFrames_;
    return true;
}

void Trajectory::read_file_step0(std::string filename) {
    // read all into cache
    ifstream ifs(filename);
    if(not ifs)
        ERROR("Can't open ["+filename+"] to read.");
    int lineno = 0;
    string line;
    this->trajFile.filename = filename;
    this->trajFile.lines.clear();
    output("Reading ["+filename+"]...");
    auto start = chrono::system_clock::now();
    while(getline(ifs,line)){
        trajFile.lines.push_back(line);
        if(++lineno > MY_FILE_LINES_UPPER_BOUND)
            ERROR("\nFile ["+filename+"] too large (>"+
                  to_string(MY_FILE_LINES_UPPER_BOUND)+") lines, reading aborted!");
    }
    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
    ostringstream oss;
    oss<<"Done in "<<duration.count()/1000.0<<" seconds."<<endl;
    output(oss.str());
}
void Trajectory::read_preparation_step1_find_frames_in_file(){
    int lineno = 0;
    trajFile.timesteps.clear();
    trajFile.startlines.clear();
    try {
        while (lineno < trajFile.lines.size()) {
            trajFile.startlines.push_back(lineno);
            my_assert_true(StringStartsWith(trajFile.lines[lineno++], "ITEM: TIMESTEP"));
            int timeStep = stoi(trajFile.lines[lineno++]);
            trajFile.timesteps.push_back(timeStep);
            my_assert_true(StringStartsWith(trajFile.lines[lineno++], "ITEM: NUMBER OF ATOMS"));
            int nAtoms = stoi(trajFile.lines[lineno++]);
            // no assertion for nAtoms, for now
            my_assert_true(StringStartsWith(trajFile.lines[lineno++], "ITEM: BOX BOUNDS"));
            lineno+=3;
            my_assert_true(StringStartsWith(trajFile.lines[lineno++], "ITEM: ATOMS"));
            lineno+=nAtoms;
        }
    }
    catch(exception e){
        ERROR("\nWhile reading line "+to_string(lineno+1)+" of ["+trajFile.filename+
              "]:\n\""+trajFile.lines[lineno]+"\"");
    }
}
void Trajectory::read_preparation_step2_find_actually_read_frames(int maxFrames, set<int> &certainFrames){
    // the [certainFrames] parameter has higher priority
    if(not certainFrames.empty()){
        vector<int> temp1, temp2;
        for(int i=0;i<trajFile.timesteps.size();i++){
            if(certainFrames.contains( trajFile.timesteps[i]) ){
                temp1.push_back(trajFile.timesteps[i]);
                temp2.push_back(trajFile.startlines[i]);
            }
        }
        trajFile.timesteps = temp1;
        trajFile.startlines = temp2;
    }
    if(trajFile.timesteps.size()>maxFrames){
        auto &v = trajFile.timesteps;
        v.erase(v.begin()+maxFrames,v.end());
        auto &w = trajFile.startlines;
        w.erase(w.begin()+maxFrames,w.end());
    }
}
void Trajectory::read_preparation_step3_remove_duplication(bool removeDup) {
    if(not removeDup)
        return;
    // Remove duplicate frames as requested. Remove from back to front
    // Assume that duplication can occur only in different files.
    if(removeDup){
        int oldNFrames = nFrames_;
        my_assert_true(oldNFrames==ts_.size());
        int deletedFrames = 0;
        for(int i=oldNFrames-1;i>=0;i--){
            bool dup = false;
            for(int j=0;j<trajFile.timesteps.size();j++){
                if(ts_[i] == trajFile.timesteps[j]){
                    dup = true;
                    break;
                }
            }
            ostringstream oss;
            oss<<"Data for timestep "<<ts_[i]<<" appears again in ["<<trajFile.filename<<"]. Old data for this timestep is removed.";
            WARNING(oss.str());
            this->DiscardFrame(i);
            deletedFrames++;
            my_assert_true(nFrames_ + deletedFrames==oldNFrames);
        }
    }
}

void Trajectory::read_one_frame(int iFrame,int iFrameInTrajFile){
    int startline = trajFile.startlines[iFrameInTrajFile];
    int lineno = startline;
    int nAtoms = stoi(trajFile.lines[lineno+3]);

    // Write systemwise info
    ts_[iFrame] = trajFile.timesteps[iFrameInTrajFile];

    // Allocate per-atom memory
    while(not StringRegexMatch(trajFile.lines[++lineno], "ITEM: ATOMS"));
    KeywordsColumnPos kcp;
    kcp.FindColumnPos(trajFile.lines[lineno]);

    if(kcp.id==-1)
        ERROR("\nTrajectory file ["+trajFile.filename+"] missing key info: id");
    if(kcp.x==-1 or kcp.y==-1 or kcp.z==-1)
        ERROR("\nTrajectory file ["+trajFile.filename+"] missing key info: coordinates");

    this->createMemoryForFrame(iFrame,kcp);
    // Read every atom,
    for(int i=0;i<nAtoms;i++){
        auto parts = StringSplit(trajFile.lines[++lineno]);
        int index = stoi(parts[kcp.id])-1;
        this->x_[iFrame][index][0] = stof(parts[kcp.x]);
        this->x_[iFrame][index][1] = stof(parts[kcp.y]);
        this->x_[iFrame][index][2] = stof(parts[kcp.z]);
        // The following info are not always present
        if(kcp.vx!=-1 or kcp.vy!=-1 or kcp.vz!=-1) {
            // if any one of vx, vy, or vz is present, the entire space for v is created.
            this->v_[iFrame][index][0] = stof(parts[kcp.vx]);
            this->v_[iFrame][index][1] = stof(parts[kcp.vy]);
            this->v_[iFrame][index][2] = stof(parts[kcp.vz]);
        }
        if(kcp.fx!=-1 or kcp.fy!=-1 or kcp.fz!=-1) {
            // if any one of vx, vy, or vz is present, the entire space for v is created.
            this->f_[iFrame][index][0] = stof(parts[kcp.fx]);
            this->f_[iFrame][index][1] = stof(parts[kcp.fy]);
            this->f_[iFrame][index][2] = stof(parts[kcp.fz]);
        }
        if(kcp.ix!=-1 or kcp.iy!=-1 or kcp.iz!=-1) {
            // if any one of vx, vy, or vz is present, the entire space for v is created.
            this->i_[iFrame][index][0] = stof(parts[kcp.ix]);
            this->i_[iFrame][index][1] = stof(parts[kcp.iy]);
            this->i_[iFrame][index][2] = stof(parts[kcp.iz]);
        }
    }

}

int Trajectory::Read(string filename, int max_workers, int maxFrames, bool removeDup, set<int> certainFrames){
    read_file_step0(filename);
    read_preparation_step1_find_frames_in_file();
    read_preparation_step2_find_actually_read_frames(maxFrames,certainFrames);
    read_preparation_step3_remove_duplication(removeDup);

    // Now resize all vectors: Note that the expanded vectors are not initialized until
    // the moment each frame is read in.
    int oldNFrames = nFrames_;
    nFrames_ += trajFile.timesteps.size();
    OPERATION_FOR_ALL_SYSTEM_VECTORS(vector_resize,nFrames_);
    OPERATION_FOR_ALL_PERATOM_VECTORS(vector_resize,nFrames_);

    for(int i=oldNFrames;i<nFrames_;i++){
        read_one_frame(oldNFrames+i, i);
    }
    return 0;
}


void Trajectory::_thread_main_(int iThread,int iFrameStart, int iFrameEnd,promise<bool> ret) {

    uniform_int_distribution<int> uid(10,50);

    for(int iFrame=iFrameStart;iFrame<iFrameEnd;iFrame++){
        /* Each Thread work independently, no need to add lock? */
        KeywordsColumnPos kcp;
        kcp.x = 0;
        createMemoryForFrame(iFrame,kcp);
        for(int iAtom=0;iAtom<nAtoms_;iAtom++){
            x_[iFrame][iAtom] = {double(iFrame),double(iAtom),double(iFrame)*iAtom};
        }
    }
    this_thread::sleep_for(chrono::milliseconds(uid(e))); // sleep for a random time
    ret.set_value(true);
}


void Trajectory::_testMultiThread() {
    nAtoms_ = 114514;
    nFrames_ = 1231;
    int nThreads = 8;

    OPERATION_FOR_ALL_SYSTEM_VECTORS(vector_resize,nFrames_);
    OPERATION_FOR_ALL_PERATOM_VECTORS(vector_resize,nFrames_);

    nThreads = min(8,nFrames_);
    int framesEachThread = nFrames_/nThreads;

    vector<future<bool>> rets;
    for(int iThread=0;iThread<nThreads;iThread++){
        int start = iThread*framesEachThread;
        int end = (iThread==nThreads-1? nFrames_ : (iThread+1)*framesEachThread);
        promise<bool> ret_state;
        rets.push_back(ret_state.get_future());
        thread th(&Trajectory::_thread_main_,this,iThread,start,end,move(ret_state));
        th.detach();
        //(iThread,start,end);
        cout<<"Thread "<<iThread<<" started"<<endl;

    }

    for(int iThread=0;iThread<nThreads;iThread++){
        bool state = rets[iThread].get();
        cout<<"Thread "<<iThread<<" ended : "<<state<<endl;
    }

    // check consistency

    for(int iFrame=0;iFrame<nFrames_;iFrame++){
        for(int iAtom=0;iAtom<nAtoms_;iAtom++){
            MY_ASSERT_DOUBLE_EQUAL(x_[iFrame][iAtom][0], iFrame);
            MY_ASSERT_DOUBLE_EQUAL(x_[iFrame][iAtom][1], iAtom);
            MY_ASSERT_DOUBLE_EQUAL(x_[iFrame][iAtom][2], double(iFrame)*iAtom);
        }
        ProgressBar(double(iFrame)/nFrames_);
    }
    ProgressBar(1.0);
}
