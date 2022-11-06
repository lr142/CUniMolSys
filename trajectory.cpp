#include "trajectory.h"
#include "opener.h"
#include "molecularmanipulator.h"
#include <future>
#include <fstream>
#include <string>
#include <sstream>
#include <memory>
#include <algorithm>
using namespace std;

// Auxiliary functions used in this module

template<class T> void ptr_create_space(T* &ptr, int size, bool nullptr_condition){
    if(nullptr_condition)
        ptr = nullptr;
    else{
        try{
            ptr = new T[size];
        }
        catch(exception){
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

int TrajFile::NFrames(){
    ERROR_IF_FALSE(timesteps.size() == lines.size())
    ERROR_IF_FALSE(timesteps.size() == nAtoms.size())
    return timesteps.size();
}
void TrajFile::Clear() {
    timesteps.clear();
    lines.clear();
    nAtoms.clear();
    filename = "";
}
void TrajFile::AddFrame(int nAtoms,long long ts){
    timesteps.push_back(ts);
    this->nAtoms.push_back(nAtoms);
    lines.push_back(make_shared<vector<string>>());
    lines[lines.size()-1]->reserve(nAtoms+1); // one extra line for the "ITEM: ATOMS"
}

void KeywordsColumnPos::FindColumnPos(string line) {
    // a line like this: "ITEM: ATOMS id mol type x y z vx vy vz"
    auto parts = StringSplit(line,' ');
    parts.erase(parts.begin(),parts.begin()+2);
    map<string,int> k2c; // key to column pos map
    for(unsigned int i=0;i<parts.size();i++){
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
TrajectoryFrame::TrajectoryFrame(int nAtoms){
    nAtoms_ = nAtoms;
    ts_ = 0;
    s_ = nullptr;
    x_ = v_ = f_ =  nullptr;
    i_ = nullptr;
}
TrajectoryFrame::~TrajectoryFrame(){
    destroyMemory();
    //output("Frame with ts = "+to_string(this->ts_)+" destoryed.");
}
void TrajectoryFrame::Read(TrajFile &trajFile, int iFrameInTrajFile) {
    // Write systemwise info
    this->ts_ = trajFile.TimeStep(iFrameInTrajFile);
    // double-checking.
    ERROR_IF_FALSE(this->nAtoms_ == trajFile.NAtoms(iFrameInTrajFile));

    // Allocate per-atom memory
    auto lines = trajFile[iFrameInTrajFile];
    KeywordsColumnPos kcp;
    ERROR_IF_FALSE(StringRegexMatch(lines[0], "ITEM: ATOMS"));
    kcp.FindColumnPos(lines[0]);

    if(kcp.id==-1)
        ERROR("\nTrajectory file ["+trajFile.filename+"] missing key info: id");
    if(kcp.xu==-1 or kcp.yu==-1 or kcp.zu==-1)
        ERROR("\nTrajectory file ["+trajFile.filename+"] missing key info: coordinates");

    /* destroy at first in case of "re-reading" a frame at the same pos, save memory. If the memory has not
    yet created, nothing will happen */
    this->destroyMemory();
    this->createMemory(kcp);

    // Read every atom, Starting from the second line, each line represents an atom.
    int i;
    try {
        for (i = 0; i < nAtoms_; i++) {
            auto parts = StringSplit(lines[i + 1]);
            int serial = stoi(parts[kcp.id]);
            // serial and xyz must always be present
            s_[i] = serial;
            x_[i][0] = stof(parts[kcp.xu]);
            x_[i][1] = stof(parts[kcp.yu]);
            x_[i][2] = stof(parts[kcp.zu]);
            // The following info are not always present
            if (kcp.vx != -1 or kcp.vy != -1 or kcp.vz != -1) {
                // if any one of vx, vy, or vz is present, the entire space for v is created.
                // Same rules apply also for fx~z, ix~z
                v_[i][0] = stof(parts[kcp.vx]);
                v_[i][1] = stof(parts[kcp.vy]);
                v_[i][2] = stof(parts[kcp.vz]);
            }
            if (kcp.fx != -1 or kcp.fy != -1 or kcp.fz != -1) {
                f_[i][0] = stof(parts[kcp.fx]);
                f_[i][1] = stof(parts[kcp.fy]);
                f_[i][2] = stof(parts[kcp.fz]);
            }
            if (kcp.ix != -1 or kcp.iy != -1 or kcp.iz != -1) {
                i_[i][0] = stoi(parts[kcp.ix]);
                i_[i][1] = stoi(parts[kcp.iy]);
                i_[i][2] = stoi(parts[kcp.iz]);
            }
        }
    }
    catch(exception){
        ERROR("At this line:\n"+lines[i + 1]);
    }
    sort_atoms();
    // construct the serial_to_index map
    this->s_to_index.clear();
    for(int i=0;i<nAtoms_;i++){
        s_to_index[s_[i]] = i;
    }
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
    sort(vm.begin(),vm.end(),[](auto &item1, auto &item2){ return item1.first < item2.first;});
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

//Trajectory
Trajectory::Trajectory(MolecularSystem &ms):ms_(ms){}
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
int Trajectory::NFrames() {
    return frames_.size();
}

bool read_one_frame_from_file_(ifstream &ifs, int &lineno,TrajFile& trajFile,set<long long> &certainFrames){
    string line;
    int nAtoms;
    long long ts;
    // Read the header
    // Note: lineno is incremented AFTER the line is read and parsed.
    if(not getline(ifs,line) or not StringRegexMatch(line, "ITEM: TIMESTEP"))
        return false;
    lineno++;
    getline(ifs,line);
    ts = stoll(StringSplit(line)[0]);
    lineno++;
    getline(ifs,line);
    if(not StringRegexMatch(line,"ITEM: NUMBER OF ATOMS"))
        throw exception();
    lineno++;
    getline(ifs,line);
    nAtoms = stoi(StringSplit(line)[0]);
    lineno++;
    for(int i=0;i<4;i++){
        lineno++;
        getline(ifs,line);
    }
    if( certainFrames.empty() or certainFrames.contains(ts)) {
        // Create a frame
        trajFile.AddFrame(nAtoms, ts);
        // there are nAtoms+1 lines to add into the last frame of trajFile
        for (int i = 0; i < nAtoms + 1; i++) {
            getline(ifs, line);
            trajFile[trajFile.NFrames() - 1].push_back(line);
        }
    }else{
        // read without storing anything
        for (int i = 0; i < nAtoms + 1; i++) {
            getline(ifs, line);
        }
    }
    lineno += (nAtoms + 1);
    return true;
}

void Trajectory::read_preparation_step1_find_frames_in_file(string filename,int maxFrames, set<long long> &certainFrames){
    // read all into cache
    ifstream ifs(filename);
    if(not ifs)
        ERROR("Can't open ["+filename+"] to read.");

    string line;
    trajFile.Clear();
    trajFile.filename = filename;

    output("Reading ["+filename+"] into memory ...");
    auto start = chrono::system_clock::now();

    int lineno = 0;
    try {
        while(read_one_frame_from_file_(ifs,lineno,trajFile,certainFrames)){
            if(trajFile.NFrames() >= maxFrames)
                break;
        }
    }
    catch(exception){
        ERROR("\nWhile reading line "+to_string(lineno+1)+" of ["+trajFile.filename+"]");
    }

    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
    ostringstream oss;
    oss<<"Done in "<<duration.count()/1000.0<<" seconds. File contains "<<trajFile.NFrames()<<" frames.";
    output(oss.str());

    // double check
    for(int i=0;i<trajFile.NFrames();i++){
        if( (int)trajFile[i].size() != trajFile.NAtoms(i)+1 or
            not StringRegexMatch(trajFile[i][0],"ITEM: ATOMS"))
            ERROR("Trajectory file reading error for frame ."+ to_string(i));
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
            for(int j=0;j<trajFile.NFrames();j++){
                if(frames_[i]->ts_ == trajFile.TimeStep(j)){
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
            ERROR_IF_FALSE(NFrames()+deletedFrames==oldNFrames);
        }
    }
}
void Trajectory::read_initialize_step4_initialize_frame_vector(){
    // Now resize all vectors: Note that the expanded vectors are not initialized until
    // the moment each frame is read in.
    int oldNFrames = NFrames();
    int newNFrames = oldNFrames + trajFile.NFrames();
    frames_.resize(newNFrames);
    for(int i=0;i<trajFile.NFrames();i++){
        frames_[oldNFrames+i] =
                make_shared<TrajectoryFrame>(trajFile.NAtoms(i));
    }
}
void Trajectory::read_final_step5_in_parallel(int oldNFrames, int max_workers){
    // Read all frames, in parallel, with timing...
    output("Processing "+to_string(trajFile.NFrames())+" new frames..");
    auto start = chrono::system_clock::now();
    if(max_workers<1) // automatically determine threads number;
        max_workers = thread::hardware_concurrency();
    vector<thread> th_;
    for(int iThread=0;iThread<max_workers;iThread++){
        th_.push_back(thread(&Trajectory::__read_frame_thread_main__, this, iThread, max_workers, oldNFrames));
    }
    for(int iThread=0;iThread<max_workers;iThread++){
        th_[iThread].join();
    }
    ProgressBar(1.0);
    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
    output("Done in "+ to_string(duration.count()/1000.0)+" seconds.");
}


int Trajectory::Read(string filename, int max_workers, int maxFrames, bool removeDup, set<long long> certainFrames){


    read_preparation_step1_find_frames_in_file(filename,maxFrames,certainFrames);
    // read_preparation_step2_find_actually_read_frames(maxFrames,certainFrames);
    read_preparation_step3_remove_duplication(removeDup);
    int oldNFrames = NFrames(); // initialize oldNFrame after step2 because frames may be deleted in step2.
    read_initialize_step4_initialize_frame_vector();
    bool testing = false;
    if(testing){ // for testing scalability
        output(" 1 worker");
        read_final_step5_in_parallel(oldNFrames, 1);
        output(" 2 workers");
        read_final_step5_in_parallel(oldNFrames, 2);
        output(" 4 workers");
        read_final_step5_in_parallel(oldNFrames, 4);
        output(" 8 workers");
        read_final_step5_in_parallel(oldNFrames, 8);
    }else{
        read_final_step5_in_parallel(oldNFrames, max_workers);
    }
    return NFrames()-oldNFrames;
}

void Trajectory::__read_frame_thread_main__(int iThread, int NThreads, int oldNFrames) {
    int nTasks = NFrames() - oldNFrames;
    auto pair = TaskDistribution(iThread,NThreads,nTasks);
    int start = pair.first;
    int end = pair.second;

    for(int i=start; i<end; i++){
        //cout<<"Thread "<<iThread<<" Prasing Frame "<<localIFrame<<"/"<<Nframes<<"..."<<endl;
        frames_[i+oldNFrames]->Read(trajFile,i);
        //cout<<"Thread "<<iThread<<" Parsed Frame "<<localIFrame<<"/"<<Nframes<<endl;
        if(iThread==0)
            ProgressBar(1.0*i/(end-start));
    }
}


shared_ptr<MolecularSystem> Trajectory::UpdateCoordsAtFrame(int iFrame, const set<int>& update_indexes, const set<int>& report_indexes){
    /* This function should be thread-safe. It only reads member data, no writes on shared member data */
    ostringstream err_msg;
    if(iFrame<0 or iFrame>NFrames()){
        err_msg << "iFrame = " << iFrame << " out of range (1~" << NFrames() << ")";
        ERROR(err_msg.str());
    }
    TrajectoryFrame& f = (*this)[iFrame];

    MolecularSystem copy = ms_.DeepCopy();
    MolecularSystemAccessor msa_copy(copy);

    vector<bool> mask(msa_copy.AtomsCount(),false);
    for(int i=0;i<msa_copy.AtomsCount();i++){
        if(report_indexes.empty() or report_indexes.contains(i)){
            mask[i] = true;
        }
    }
    for(int i=0;i<msa_copy.AtomsCount();i++){
        if(mask[i] and (update_indexes.empty() or update_indexes.contains(i))){
            int serial = i+1;
            int index_in_traj = f.SerialToIndex(serial);
            if(index_in_traj == -1)
                continue;
            msa_copy.AtomByGlobalIndex(i).xyz = f.X()[index_in_traj];
        }
    }
    MolSysSubsystemByMask(copy,mask);
    shared_ptr<MolecularSystem> pMS = make_shared<MolecularSystem>(copy);
    pMS->SetName(ms_.GetName()+" "+to_string(this->frames_[iFrame]->ts_));
    return pMS;
}

void Trajectory::__show_trajectory_thread_main__(int iThread,int NThreads,const set<int>& update_indexes,
                                                 const set<int>& report_indexes, vector<shared_ptr<MolecularSystem>> &sys){
    int nTasks = NFrames();
    auto pair = TaskDistribution(iThread,NThreads,nTasks);
    for(int i=pair.first;i<pair.second;i++){
        sys[i] = UpdateCoordsAtFrame(i,update_indexes,report_indexes);
        MolSysReduceToSingleMolecule(*sys[i]);
        (*sys[i])[0].name = to_string(frames_[i]->ts_);
        if(iThread==0)
            ProgressBar(1.0*i/(pair.second-pair.first));
    }
    if(iThread==0)
        ProgressBar(1.0);
}

void Trajectory::ShowTrajectory(string filename, const set<int>& update_indexes, const set<int>& report_indexes, int max_workers) {
    vector<shared_ptr<MolecularSystem>> molsys_frames(NFrames());
    if(max_workers < 1)
        max_workers = thread::hardware_concurrency();
    max_workers = min(max_workers,NFrames());

    output("Generate trajectory for "+ to_string(NFrames())+" frames...");
    auto start = chrono::system_clock::now();

    vector<thread> th_;
    for(int i=0;i<max_workers;i++){
        th_.push_back(thread(&Trajectory::__show_trajectory_thread_main__,this,i,max_workers,
                             ref(update_indexes),ref(report_indexes),ref(molsys_frames)));
    }
    MolecularSystem entire;
    for(int i=0;i<(int)th_.size();i++){ // must join in order
        th_[i].join();
        auto pair = TaskDistribution(i,max_workers,NFrames());
        for(int j=pair.first;j<pair.second;j++)
            entire.AddMolecule((*molsys_frames[j])[0]);
    }
    if(filename != "")
        QuickSave(entire,filename);

    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
    output("Done in "+ to_string(duration.count()/1000.0) + " seconds.");
}
