#include "trajectory.h"
Trajectory::Trajectory(MolecularSystem &ms):ms_(ms){
    nAtoms_ = ms.AtomsCount();
    nFrames_ = 0;
}
Trajectory::~Trajectory(){
    destroyAll();
}
void Trajectory::createMemoryForFrame(int iFrame) {
#define CREATE_ONE_OBJ(OBJ,T) \
    if((OBJ).size()!=0){ (OBJ)[iFrame] = new XYZ_T_<T>[nAtoms_]; }
    CREATE_ONE_OBJ(x_,double)
    CREATE_ONE_OBJ(v_,double)
    CREATE_ONE_OBJ(f_,double)
    CREATE_ONE_OBJ(image_flags_,int)
}
void Trajectory::destroyMemoryForFrame(int iFrame){
#define DESTROY_ONE_OBJ(OBJ) \
    if((OBJ).size()!=0){ delete [] (OBJ)[iFrame]; }
    DESTROY_ONE_OBJ(x_)
    DESTROY_ONE_OBJ(v_)
    DESTROY_ONE_OBJ(f_)
    DESTROY_ONE_OBJ(image_flags_)
}
void Trajectory::destroyAll() {
    for(int i=0;i<nFrames_;i++){
        destroyMemoryForFrame(i);
    }
    x_.clear();
    v_.clear();
    f_.clear();
    image_flags_.clear();
}
