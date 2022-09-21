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

template<class T> void ptr_create_space(T* &ptr, int size, bool nullptr_condition){
    if(nullptr_condition)
        ptr = nullptr;
    else
        ptr = new T[size];
}
template <class T> void ptr_destroy_space(T* &ptr){
    if(ptr != nullptr){
        delete [] ptr;
        ptr = nullptr;
    }
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
    // Special rules: ideally we want to read xu,yu,zu.
    // However, if xu (unwrapped) not reported but x (wrapped), xs(scaled), xsu(scaled unwrapped) is reported
    // Use that in the line of succession xu <-- x<-- xs <-- xsu (a left value has higher priority than the right,
    // use the right value only if the left value is unavailable. )
#define REPLACE_IF(VAR1, VAR2) if (VAR1==-1 and VAR2!=-1){VAR1=VAR2;}
    // The replacement chain must work in reverse.
    REPLACE_IF(xs,xsu)  REPLACE_IF(x,xs)  REPLACE_IF(xu,x)
    REPLACE_IF(ys,ysu)  REPLACE_IF(y,ys)  REPLACE_IF(yu,y)
    REPLACE_IF(zs,zsu)  REPLACE_IF(z,zs)  REPLACE_IF(zu,z)
#undef REPLACE_IF
}

// TrajectoryFrame
TrajectoryFrame::TrajectoryFrame() {
    nAtoms_ = ts_ = 0;
    s_ = nullptr;
    x_ = v_ = f_ =  nullptr;
    i_ = nullptr;
}
TrajectoryFrame::~TrajectoryFrame(){
    destroyMemory();
}
void TrajectoryFrame::Read(TrajFile &trajFile,int iFrameInTrajFile){
    int startline = trajFile.startlines[iFrameInTrajFile];
    int lineno = startline;

    // Write systemwise info
    ts_ = trajFile.timesteps[iFrameInTrajFile];
    nAtoms_ = stoi(trajFile.lines[lineno+3]);

    // Allocate per-atom memory
    JumpToLine(trajFile.lines,"ITEM: ATOMS",lineno,startline,startline+20);
    KeywordsColumnPos kcp;
    kcp.FindColumnPos(trajFile.lines[lineno]);

    if(kcp.id==-1)
        ERROR("\nTrajectory file ["+trajFile.filename+"] missing key info: id");
    if(kcp.xu==-1 or kcp.yu==-1 or kcp.zu==-1)
        ERROR("\nTrajectory file ["+trajFile.filename+"] missing key info: coordinates");

    this->createMemory(kcp);
    // Read every atom,
    for(int i=0;i<nAtoms_;i++){
        auto parts = StringSplit(trajFile.lines[++lineno]);
        int serial = stoi(parts[kcp.id]);
        // serial and xyz must always be present
        s_[i] = serial;
        x_[i][0] = stof(parts[kcp.xu]);
        x_[i][1] = stof(parts[kcp.yu]);
        x_[i][2] = stof(parts[kcp.zu]);
        // The following info are not always present
        if(kcp.vx!=-1 or kcp.vy!=-1 or kcp.vz!=-1) {
            // if any one of vx, vy, or vz is present, the entire space for v is created.
            v_[i][0] = stof(parts[kcp.vx]);
            v_[i][1] = stof(parts[kcp.vy]);
            v_[i][2] = stof(parts[kcp.vz]);
        }
        if(kcp.fx!=-1 or kcp.fy!=-1 or kcp.fz!=-1) {
            // if any one of vx, vy, or vz is present, the entire space for v is created.
            f_[i][0] = stof(parts[kcp.fx]);
            f_[i][1] = stof(parts[kcp.fy]);
            f_[i][2] = stof(parts[kcp.fz]);
        }
        if(kcp.ix!=-1 or kcp.iy!=-1 or kcp.iz!=-1) {
            // if any one of vx, vy, or vz is present, the entire space for v is created.
            i_[i][0] = stof(parts[kcp.ix]);
            i_[i][1] = stof(parts[kcp.iy]);
            i_[i][2] = stof(parts[kcp.iz]);
        }
    }
    sort_atoms();
}
void TrajectoryFrame::createMemory(KeywordsColumnPos &kcp) {
    bool cond;
    if(nAtoms_ == 0)
        ERROR("nAtoms_ == 0, Can't createMemory.");
    // Unfortunately, we can't write a macro to make this function simpler.
    ptr_create_space(s_,nAtoms_,false); // always create this

    ptr_create_space(x_,nAtoms_,false); // always create this

    cond = kcp.vx<0 and kcp.vy<0 and kcp.vz<0;
    ptr_create_space(v_,nAtoms_,cond);

    cond = kcp.fx<0 and kcp.fy<0 and kcp.fz<0;
    ptr_create_space(f_,nAtoms_,cond);

    cond = kcp.ix<0 and kcp.iy<0 and kcp.iz<0;
    ptr_create_space(i_,nAtoms_,cond);
}
void TrajectoryFrame::destroyMemory(){
    ptr_destroy_space(s_);
    ptr_destroy_space(x_);
    ptr_destroy_space(v_);
    ptr_destroy_space(f_);
    ptr_destroy_space(i_);
}
template <class T> void sort_array(T* array, int* key_array, int length){
    // will keep key_array unchanged.
    if(array== nullptr)
        return;
    vector<pair<int,T>> vm(length);
    for(int i=0;i<length;i++){
        vm[i] = make_pair(key_array[i],array[i]);
    }
    ::sort(vm.begin(),vm.end(),[](auto &item1, auto &item2){ return item1.first < item2.first;});
    for(int i=0;i<length;i++){
        array[i] = vm[i].second;
    }
}
void TrajectoryFrame::sort_atoms(){
    if(nAtoms_ == 0)
        return;
    // sort_atoms according to serials s_[]
    sort_array(x_,s_,nAtoms_);
    sort_array(v_,s_,nAtoms_);
    sort_array(f_,s_,nAtoms_);
    sort_array(i_,s_,nAtoms_);
    // must sort_atoms s_ itself at last.
    ::sort(s_, s_ + nAtoms_);
    /* Debug, show atoms */
//    for(int i=0;i<nAtoms_;i++){
//        cout<<i<<" "<<s_[i]<<" "<<x_[i]<<" "<<v_[i]<<endl;
//    }
}

//Trajector
Trajectory::Trajectory(MolecularSystem &ms):ms_(ms),e(0){
}
Trajectory::~Trajectory() {
    Clear();
}
bool Trajectory::DiscardFrame(int iFrame){
    if(iFrame<0 or iFrame>=NFrames())
        return false;
    frames_.erase(frames_.begin()+iFrame,frames_.begin()+iFrame+1);
    return true;
}
void Trajectory::Clear() {
    frames_.clear();
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
        int oldNFrames = NFrames();
        int deletedFrames = 0;
        for(int i=oldNFrames-1;i>=0;i--){
            bool dup = false;
            for(int j=0;j<trajFile.timesteps.size();j++){
                if(frames_[i].ts_ == trajFile.timesteps[j]){
                    dup = true;
                    break;
                }
            }
            ostringstream oss;
            oss<<"Data for timestep "<<frames_[i].ts_<<" appears again in ["<<trajFile.filename<<"]. Old data for this old frame will be removed.";
            WARNING(oss.str());
            this->DiscardFrame(i);
            deletedFrames++;
            my_assert_true(NFrames()+deletedFrames==oldNFrames);
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
    int oldNFrames = NFrames();
    int newNFrames = oldNFrames + trajFile.timesteps.size();

    frames_.resize(newNFrames);

    for(int i=0;i<trajFile.timesteps.size();i++){
        frames_[oldNFrames+i].Read(trajFile,i);
    }
    return 0;
}


void Trajectory::_thread_main_(int iThread,int iFrameStart, int iFrameEnd, int NAtoms) {
    uniform_int_distribution<int> uid(10,50);
    for(int iFrame=iFrameStart;iFrame<iFrameEnd;iFrame++){
        /* Each Thread work independently, no need to add lock? */
        KeywordsColumnPos kcp;
        kcp.x = 0;
        auto &frame = frames_[iFrame];
        frame.nAtoms_ = NAtoms;
        frame.createMemory(kcp);
        for(int iAtom=0;iAtom<NAtoms;iAtom++){
            frame.x_[iAtom] = {double(iFrame),double(iAtom),double(iFrame)*iAtom};
        }
    }
    this_thread::sleep_for(chrono::milliseconds(uid(e))); // sleep for a random time
}


bool double_equal(double a,double b){
    return fabs(a-b)<MY_SMALL;
}
void Trajectory::_testMultiThread() {
    int NAtoms = 114514;
    int NFrames = 964;
    int nThreads = 8;

    frames_.resize(NFrames);

    nThreads = min(8,NFrames);
    int framesEachThread = NFrames/nThreads;

    vector<thread> th_;
    for(int iThread=0;iThread<nThreads;iThread++){
        int start = iThread*framesEachThread;
        int end = (iThread==nThreads-1? NFrames : (iThread+1)*framesEachThread);
        th_.push_back(thread(&Trajectory::_thread_main_,this,iThread,start,end,NAtoms));
        //th_[th_.size()-1].detach();
        //(iThread,start,end);
        cout<<"Thread "<<iThread<<" started"<<endl;
    }

    for(int iThread=0;iThread<nThreads;iThread++){
        th_[iThread].join();
        cout<<"Thread "<<iThread<<" ended : "<<endl;
    }

    // check consistency

    for(int iFrame=0;iFrame<NFrames;iFrame++){
        for(int iAtom=0;iAtom<NAtoms;iAtom++){
#define MY_ASSERT_DOUBLE_EQUAL(a,b) if(not double_equal((a),(b))){cout<<iFrame<<" "<<iAtom<<" "<<(a)<<" "<<(b)<<endl;}
            MY_ASSERT_DOUBLE_EQUAL(frames_[iFrame].x_[iAtom][0], iFrame);
            MY_ASSERT_DOUBLE_EQUAL(frames_[iFrame].x_[iAtom][1], iAtom);
            MY_ASSERT_DOUBLE_EQUAL(frames_[iFrame].x_[iAtom][2], double(iFrame)*iAtom);
#undef MY_ASSERT_DOUBLE_EQUAL
        }
        if(nThreads==1)
            ProgressBar(double(iFrame)/iFrame);
    }
    ProgressBar(1.0);
}
