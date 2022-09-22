#include "trajectory.h"
#include "opener.h"
#include "molecularmanipulator.h"
#include <future>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
using namespace std;

// Auxiliary functions used in this module
void my_assert_true(bool statement){
    if(not (statement))
        throw exception();
}

template<class T> void ptr_create_space(T* &ptr, int size, bool nullptr_condition){
    if(nullptr_condition)
        ptr = nullptr;
    else{
        try{
            ptr = new T[size];
        }
        catch(exception e){
            ERROR("Memory Allocation Failed.");
        }
    }
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
    //output("Frame with ts = "+to_string(this->ts_)+" destoryed.");
}
void TrajectoryFrame::Read(TrajFile &trajFile, int iFrameInTrajFile, bool createMemory) {
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

    if(createMemory)
        this->createMemory(kcp);

    // Read every atom,
    for(int i=0;i<nAtoms_;i++){
        auto parts = StringSplit(trajFile.lines[++lineno]);
        int serial = stoi(parts[kcp.id]);
        // serial and xyz must always be present
        s_[i] = serial;
        if(s_[i] < 1){
            ostringstream oss;
            oss<<"\nTrajectory file ["<<trajFile.filename<<" @ line "<<lineno+1<<"\n"
            <<trajFile.lines[lineno]<<"\n"<<"atom id can't be < 1";
            ERROR(oss.str());
        }
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
            i_[i][0] = stoi(parts[kcp.ix]);
            i_[i][1] = stoi(parts[kcp.iy]);
            i_[i][2] = stoi(parts[kcp.iz]);
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
    output("Reading ["+filename+"] into memory ...");
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
    oss<<"Done in "<<duration.count()/1000.0<<" seconds.";
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
                if(frames_[i]->ts_ == trajFile.timesteps[j]){
                    dup = true;
                    break;
                }
            }
            if(dup) {
                ostringstream oss;
                oss << "Data for timestep " << frames_[i]->ts_ << " appears again in [" << trajFile.filename
                    << "]. Old data for this frame will be removed.";
                WARNING(oss.str());
                this->DiscardFrame(i);
                deletedFrames++;
            }
            my_assert_true(NFrames()+deletedFrames==oldNFrames);
        }
    }
}
void Trajectory::read_initialize_step4_initialize_frame_vector(int &oldNFrames){
    // Now resize all vectors: Note that the expanded vectors are not initialized until
    // the moment each frame is read in.
    oldNFrames = NFrames();
    int newNFrames = oldNFrames + trajFile.timesteps.size();
    frames_.resize(newNFrames);
    /* Initialize vector and shared_ptrs. Allocate memory also at this time.
     * Note : this is to assume all frames having the same stored info. If this is not the case, you'll need to change
     * the code: move createMemory into TrajectoryFrame::Read()
     * */
    int lineno;
    int nAtoms;
    JumpToLine(trajFile.lines,"ITEM: NUMBER OF ATOMS",lineno,0,20);
    nAtoms = stoi(StringSplit(trajFile.lines[lineno+1])[0]);
    JumpToLine(trajFile.lines,"ITEM: ATOMS",lineno,0,20);
    KeywordsColumnPos kcp;
    kcp.FindColumnPos(trajFile.lines[lineno]);
    for(int i=0;i<trajFile.timesteps.size();i++){
        frames_[oldNFrames+i] = make_shared<TrajectoryFrame>();
        frames_[oldNFrames+i]->nAtoms_ = nAtoms;
//        frames_[oldNFrames+i]->createMemory(kcp);
    }
}
void Trajectory::read_final_step5_in_parallel(int oldNFrames, int max_workers){
    // Read all frames, in parallel, with timing...
    output("Processing "+to_string(trajFile.timesteps.size())+" new frames..");
    auto start = chrono::system_clock::now();
    if(max_workers<1) // automatically determine threads number;
        max_workers = thread::hardware_concurrency();
    vector<thread> th_;
    int iFrame = oldNFrames; // this is the value that's mutexed.
    for(int iThread=0;iThread<max_workers;iThread++){
        //th_.push_back(thread(&Trajectory::__read_frame_thread_main_,this,iThread,ref(iFrame),oldNFrames));
        th_.push_back(thread(&Trajectory::__read_frame_thread_main_method2_,this,iThread,max_workers,oldNFrames));
    }
    for(int iThread=0;iThread<max_workers;iThread++){
        th_[iThread].join();
    }
    ProgressBar(1.0);
    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
    output("Done in "+ to_string(duration.count()/1000.0)+" seconds.");
}


int Trajectory::Read(string filename, int max_workers, int maxFrames, bool removeDup, set<int> certainFrames){
    read_file_step0(filename);
    read_preparation_step1_find_frames_in_file();
    read_preparation_step2_find_actually_read_frames(maxFrames,certainFrames);
    read_preparation_step3_remove_duplication(removeDup);
    int oldNFrames;
    read_initialize_step4_initialize_frame_vector(oldNFrames);
    { // for testing scalability
        output(" 1 worker");
        read_final_step5_in_parallel(oldNFrames, 1);
        output(" 2 workers");
        read_final_step5_in_parallel(oldNFrames, 2);
        output(" 4 workers");
        read_final_step5_in_parallel(oldNFrames, 4);
        output(" 8 workers");
        read_final_step5_in_parallel(oldNFrames, 8);
    }

    return 0;
}
void Trajectory::__read_frame_thread_main_(int iThread,int &iFrame,int oldNFrames){
    while(true){
        int localIFrame;
        {
            mux.lock();
            //unique_lock<mutex> lock(mux);  // must add lock before read/write iFrame
            if(iFrame<NFrames()) {
                localIFrame = iFrame++;
                mux.unlock();
            }
            else{
                mux.unlock();
                break;
            }
            // it will release the lock when goes outside the bracket
        }
        // localIFrame will be defined if reached here. localIFrame-oldNFrames == iframe in the file.
        // cout<<"Thread "<<iThread<<" is reading frame "<<localIFrame<<endl;
        /* if setting the last parameter to be false, we're telling the function not to create memory */
        frames_[localIFrame]->Read(trajFile, localIFrame - oldNFrames, true);
//        if(iThread==0)
//            ProgressBar(localIFrame*1.0/NFrames());
    }
}

void Trajectory::__read_frame_thread_main_method2_(int iThread, int NThreads, int oldNFrames) {
    int Nframes = NFrames();
    int partition = Nframes%NThreads==0? Nframes/NThreads : Nframes/NThreads+1 ;
    for(int localIFrame=partition*iThread;localIFrame< min(partition*(iThread+1),Nframes);localIFrame+=1){
        frames_[localIFrame]->Read(trajFile, localIFrame - oldNFrames, true);
    }
}


shared_ptr<MolecularSystem> Trajectory::UpdateCoordsAtFrame(int iFrame, bool only_in_traj){
    ostringstream err_msg;
    if(iFrame<0 or iFrame>NFrames()){
        err_msg << "iFrame = " << iFrame << " out of range (1~" << NFrames() << ")";
        ERROR(err_msg.str());
    }
    TrajectoryFrame& f = (*this)[iFrame];

    // Find out what atoms are in the frame
    MolecularSystemAccessor msa_(ms_);
    vector<int> ms_to_traj_index_map(msa_.AtomsCount(),-1);
    for(int i=0;i<f.nAtoms_;i++){
        int index = f.s_[i] - 1; // serial to index
        /* !!!!!! If running a GCMC, possible that [index > atoms in the system]. Leave it now, just throw a exception */
        if(index<0 or index>=msa_.AtomsCount()){
            err_msg<<"Atom index read in trajectory file ("<<index<<") <0 or >= nAtoms in system ("<<msa_.AtomsCount()<<")."<<endl;
            err_msg<<"Reading GCMC trajectory is not yet implemented, sorry!";
            ERROR(err_msg.str());
        }
        // "i" is the index in trajectory, "index" is the index in original molsys;
        ms_to_traj_index_map[index] = i;
    }

    shared_ptr<MolecularSystem> pMS;
    if(only_in_traj){
        vector<bool> mask = vector<bool>(msa_.AtomsCount(),false);
        for(int i=0;i<msa_.AtomsCount();i++){
            if(ms_to_traj_index_map[i]>=0)
                mask[i] = true;
        }
        MolecularSystem copy = ms_.DeepCopy();
        MolSysSubsystemByMask(copy,mask);
        //Update coords. Atoms in pMS and in trajectory should have exactly the same atom order
        int iCounter = 0;
        for(int i=0;i<copy.MoleculesCount();i++){
            for(int j=0;j<copy[i].AtomsCount();j++){
                copy[i][j].xyz = f.x_[iCounter++];
            }
        }
        pMS = std::make_shared<MolecularSystem>(copy);
    }else{
        MolecularSystem copy = ms_.DeepCopy();
        int index = 0;
        for(int i=0;i<copy.MoleculesCount();i++){
            for(int j=0;j<copy[i].AtomsCount();j++){
                int indexInTraj = ms_to_traj_index_map[index];
                if(indexInTraj != -1)
                    copy[i][j].xyz = f.x_[indexInTraj];
                ++index;
            }
        }
    }
    return pMS;
}

void Trajectory::ShowTrajectory(std::string filename,bool showOriginal) {
    shared_ptr<MolecularSystem> pCopy;
    MolecularSystem entire;
    for(int i=0;i<NFrames();i++){
        pCopy = UpdateCoordsAtFrame(i,!showOriginal);
        MolSysReduceToSingleMolecule(*pCopy);
        entire.AddMolecule((*pCopy)[0]);
        cout<<"Frame "<<i<<endl;
    }
    QuickSave(entire,filename);
}


void Trajectory::__test_thread_main__(int iThread,int &iFrame,int nFrames, int nAtoms){
    uniform_int_distribution<int> uid(1000,4000);
    int number_of_tasks = 0;
    while(true){
        int localIFrame;
        {
            unique_lock<mutex> lock(mux);
            if(iFrame < nFrames) {
                localIFrame = iFrame;
                iFrame++;
            }
            else {
                break;
            }
        }
        cout<<"Thread "<<iThread<<" is processing task "<<iFrame<<endl;
        { // The actual work
            KeywordsColumnPos kcp;
            kcp.FindColumnPos("ITEMS: ATOM x y z id fx fy fz vx vz vy iz iy iz");
            auto &pF = frames_[localIFrame];
            pF = make_shared<TrajectoryFrame>();
            pF->nAtoms_ = nAtoms;
            pF->createMemory(kcp);
            for(int iAtom=0;iAtom<nAtoms;iAtom++){
                pF->v_[iAtom] = {XYZ_DTYPE(localIFrame),XYZ_DTYPE(iAtom),XYZ_DTYPE(localIFrame)*iAtom};
                pF->f_[iAtom] = pF->v_[iAtom] * 99.0;
                pF->x_[iAtom] = pF->f_[iAtom] * (1/99.0);
                pF->i_[iAtom] = {uid(e),uid(e),uid(e)};
            }
        }
        number_of_tasks++;
    }
    cout<<"Thread "<<iThread<<" has processed "<<number_of_tasks<<" tasks."<<endl;
}
bool double_equal(double a,double b){
    return fabs(a-b)<MY_SMALL;
}
void Trajectory::_testMultiThread() {
    int NAtoms = 114514;
    int NFrames = 1230;
    int nThreads = 8;

    frames_.resize(NFrames);

    nThreads = min(8,nThreads);

    vector<thread> th_;
    mutex mux;
    int iFrame = 0; // mutexed variable

    for(int iThread=0;iThread<nThreads;iThread++){
        th_.push_back(thread(&Trajectory::__test_thread_main__,this,iThread,ref(iFrame),NFrames,NAtoms));
        //th_[th_.size()-1].detach();
        //(iThread,start,end);
        //cout<<"Thread "<<iThread<<" started"<<endl;
    }

    for(int iThread=0;iThread<nThreads;iThread++){
        th_[iThread].join();
        cout<<"Thread "<<iThread<<" ended : "<<endl;
    }


    for(int iFrame=0;iFrame<NFrames;iFrame++){
        for(int iAtom=0;iAtom<NAtoms;iAtom++){
#define MY_ASSERT_DOUBLE_EQUAL(a,b) if(not double_equal((a),(b))){cout<<iFrame<<" "<<iAtom<<" "<<(a)<<" "<<(b)<<endl;}
            MY_ASSERT_DOUBLE_EQUAL((*this)[iFrame].x_[iAtom][0], iFrame);
            MY_ASSERT_DOUBLE_EQUAL((*this)[iFrame].x_[iAtom][1], iAtom);
            MY_ASSERT_DOUBLE_EQUAL((*this)[iFrame].x_[iAtom][2], double(iFrame)*iAtom);
#undef MY_ASSERT_DOUBLE_EQUAL
        }
        if(nThreads==1)
            ProgressBar(double(iFrame)/NFrames);
    }

    Clear();
}