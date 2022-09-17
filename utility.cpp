#include "utility.h"
#include <ctype.h>
#include <fstream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <sstream>
using namespace std;

vector<string> StringSplit(string str,char delimitator){
    string::size_type pos = 0;
    string::size_type new_pos = 0;
    vector<string> result;
    while(pos<str.size()){
        new_pos = str.find_first_of(delimitator,pos);
        if(new_pos == string::npos){
            // If the string not ends with delimitator, parse the last word
            new_pos = str.size();
        }
        if(new_pos>pos){
            // not true if two delimitator appears consecutively.
            result.push_back(str.substr(pos,new_pos-pos));
        }
        pos = new_pos + 1;
    }
    return result;
}
string StringStrip(string str){
    int beg,end;
    for(beg=0;beg<str.size();++beg){
        if(!isspace(str[beg]))
            break;
    }
    if(beg==str.size()) // all empty string
        return "";
    // The string must at least contain some non-white chars
    for(end=str.size()-1;end>=0;--end)
        if(!isspace(str[end]))
            break;
    return str.substr(beg,end-beg+1);
}
bool StringStartsWith(string str,string pattern){
    if(pattern.size() > str.size())
        return false;
    for(int i=0;i<pattern.size();++i){
        if(str[i]!=pattern[i])
            return false;
    }
    return true;
}
pair<string,string> StringTok(string str, string token){
    str = StringStrip(str);
    token = StringStrip(token);
    if(StringStartsWith(str,token)){
        auto remain = StringStrip(str.substr(token.size(),str.size()-token.size()));
        return make_pair(token,remain);
    }else{
        return make_pair("","");
    }
}
string StringToUpper(string str){
    string copy = str;
    for_each(copy.begin(),copy.end(),[](char &c){c = toupper(c);});
    return copy;
}
string StringToLower(string str){
    string copy = str;
    for_each(copy.begin(),copy.end(),[](char &c){c = tolower(c);});
    return copy;
}
string StringToCapitalized(string str){
    string res = StringToLower(str);
    res[0] = toupper(res[0]);
    return res;
}
string StringRemoveComment(string str){
    auto pos = str.find("#");
    return StringStrip(str.substr(0,pos));
}

ErrorAndOutputHandler::ErrorAndOutputHandler():onOffState_(true),pOutput_(&cout){}
void ErrorAndOutputHandler::SetOutput(std::ostream &ostream) {pOutput_ = &ostream;}
void ErrorAndOutputHandler::TurnOff() {onOffState_ = false;}
void ErrorAndOutputHandler::TurnOn() {onOffState_ = true;}

void ErrorHandler::operator()(std::string msg,const char* file,int line, bool fatal) {
    // Non-fatal messages will be turned off if the handler is in OFF state
    // But Fatal messages will be shown regardless of the on/off state
    if(!fatal and !onOffState_)
        return;
    ostringstream oss;
    oss<< (fatal? "Fatal Error [":"Warning [");
    oss<<file<<":"<<line<<"] ";
    *pOutput_<<oss.str()<<": "<<msg<<endl;
    if(fatal)
        exit(0);
}
void OutputHandler::operator()(std::string msg, bool newline) {
    if(!onOffState_)
        return;
    *pOutput_<<msg;
    if(newline)
        *pOutput_<<endl;
}

PeriodicTable::PeriodicTable() {
    // Note that all element names must be Cappitalized. Element 0 (M) is the pseudo element.
    // They must follow the order in the periodic table.
    elements_ = StringSplit("M H He "\
        "Li Be B C N O F Ne "\
        "Na Mg Al Si P S Cl Ar "\
        "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "\
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "\
        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "\
        "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo");
    string filename = DATAFILESPATH+"/atomicweights.csv";
    ifstream ifs(filename);
    if(!ifs)
        ERROR("Can't open ["+filename+"], please check the _DATAFILESPATH_ variable in the CMakeLists.txt.");
    string line;
    while(getline(ifs,line)){
        auto words = StringSplit(line,',');
        atomicWeights_[StringToCapitalized(words[0])] = stof(words[1]);
    }
//    for(auto pair:atomicWeights_){
//        cout<<pair.first<<","<<pair.second<<"|"<<endl;
//    }
    //output("PeriodicTable::PeriodicTable()");
}
PeriodicTable* PeriodicTable::getInstance_(){
    static PeriodicTable table;  // There is only one copy
    return &table;
}
string PeriodicTable::AtomicNumberToElement(int atomicNumber) {
    auto pTable = PeriodicTable::getInstance_();
    if(atomicNumber<0 || atomicNumber>pTable->elements_.size())
        return "";
    else
        return pTable->elements_[atomicNumber];
}
int PeriodicTable::ElementToAtomicNumber(string element) {
    auto pTable = PeriodicTable::getInstance_();
    element = StringToCapitalized(element);
    for(int i=0;i<pTable->elements_.size();++i){
        if(pTable->elements_[i] == element)
            return i;
    }
    return -1;
}
double PeriodicTable::AtomicWeight(string element){
    element = StringToCapitalized(element);
    auto pTable = getInstance_();
    if(pTable->atomicWeights_.find(element) == pTable->atomicWeights_.end())
        return -1.0;
    else
        return pTable->atomicWeights_[element];
}
string PeriodicTable::PossibleElementWithGivenWeight(double weight){
    auto pTable = getInstance_();
    string element_with_min_error = "";
    double min_error = 1e99;
    for(auto &item:pTable->atomicWeights_){
        double error = fabs(item.second - weight);
        if(error < min_error){
            min_error = error;
            element_with_min_error = item.first;
        }
    }
    return element_with_min_error;
}

// Mathematical Functions
XYZ::XYZ():xyz{0.0,0.0,0.0}{}
XYZ::XYZ(double x,double y,double z):xyz{x,y,z}{}
XYZ::XYZ(const double xyzs[3]):xyz{xyzs[0],xyzs[1],xyzs[2]}{}
double& XYZ::operator[](int index) {
    if(index<0 || index>2)
        throw out_of_range("XYZ[]");
    return xyz[index];
}
XYZ& XYZ::operator=(const XYZ& otherXYZ){
    for(int i=0;i<3;i++)
        xyz[i] = otherXYZ.xyz[i];
    return *this;
}
XYZ& XYZ::operator=(const double xyzs[3]){
    for(int i=0;i<3;i++)
        xyz[i] = xyzs[i];
    return *this;
}
double XYZ::NormSquared(){
    return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2];
}
double XYZ::Norm(){
    return sqrt(NormSquared());
}
XYZ& XYZ::Normalized(){
    double norm = Norm();
    if(norm < 1e-10)
        throw runtime_error("Norm too small: "+ to_string(norm));
    (*this)*=1.0/norm;
    return *this;
}
XYZ& XYZ::operator+=(XYZ v){
    xyz[0] += v[0];
    xyz[1] += v[1];
    xyz[2] += v[2];
    return *this;
}
XYZ& XYZ::operator-=(XYZ v){
    xyz[0] -= v[0];
    xyz[1] -= v[1];
    xyz[2] -= v[2];
    return *this;
}
XYZ& XYZ::operator*=(double a){
    xyz[0] *= a;
    xyz[1] *= a;
    xyz[2] *= a;
    return *this;
}
XYZ operator+(const XYZ& v1, const XYZ& v2){
    return XYZ(v1.xyz[0]+v2.xyz[0], v1.xyz[1]+v2.xyz[1], v1.xyz[2]+v2.xyz[2]);
}
XYZ operator-(const XYZ& v1, const XYZ& v2){
    return XYZ(v1.xyz[0]-v2.xyz[0], v1.xyz[1]-v2.xyz[1], v1.xyz[2]-v2.xyz[2]);
}
XYZ operator*(const XYZ& v, double a){
    return XYZ(v.xyz[0]*a, v.xyz[1]*a, v.xyz[2]*a);
}
double XYZDot(const XYZ &v1, const XYZ &v2){
    return v1.xyz[0]*v2.xyz[0] + v1.xyz[1]*v2.xyz[1] + v1.xyz[2]*v2.xyz[2];
}
XYZ XYZCross(const XYZ &v1, const XYZ &v2){
    double x =  v1.xyz[1]*v2.xyz[2] - v1.xyz[2]*v2.xyz[1];
    double y = -v1.xyz[0]*v2.xyz[2] + v1.xyz[2]*v2.xyz[0];
    double z =  v1.xyz[0]*v2.xyz[1] - v1.xyz[1]*v2.xyz[0];
    return XYZ(x,y,z);
}
bool operator==(const XYZ& v1, const XYZ& v2){
    static double tolerance = 1e-5;
    return (v1-v2).Norm() < tolerance;
}
ostream &operator<<(ostream &os, const XYZ &xyz) {
    os <<"["<<xyz.xyz[0]<<", "<<xyz.xyz[1]<<", "<<xyz.xyz[2]<<"]";
    return os;
}
void XYZRotate(XYZ* coords,int numAtoms,double degree_clockwise,XYZ axis){
    axis.Normalized();
    double rx = axis[0];
    double ry = axis[1];
    double rz = axis[2];
    // The rotational matrix was copied from a previous code. I forgot its reference...
    double phi = degree_clockwise * M_PI / 180.0;
    double c = cos(phi);
    double s = sin(phi);
    double RMatrix[3][3];
    RMatrix[0][0] = c + (1-c)*rx*rx;
    RMatrix[0][1] = (1-c)*rx*ry - rz*s;
    RMatrix[0][2] = (1-c)*rx*rz + ry*s;
    RMatrix[1][0] = (1-c)*rx*ry + rz*s;
    RMatrix[1][1] = c + (1-c)*ry*ry;
    RMatrix[1][2] = (1-c)*ry*rz - rx*s;
    RMatrix[2][0] = (1-c)*rx*rz - ry*s;
    RMatrix[2][1] = (1-c)*ry*rz + rx*s;
    RMatrix[2][2] = c+(1-c)*rz*rz;
    XYZ* newCoords = new XYZ[numAtoms];
    // newCoords = oldCoords * RMatrix, the * here is the matrix multiplication.
    for(int i=0;i<numAtoms;++i){
        for(int j=0;j<3;++j){
            newCoords[i][j] = 0.0;
            for(int k=0;k<3;++k){
                newCoords[i][j] += coords[i][k]*RMatrix[k][j];
            }
        }
    }
    for(int i=0;i<numAtoms;++i){
        coords[i] = newCoords[i];
    }
    delete [] newCoords;
}
