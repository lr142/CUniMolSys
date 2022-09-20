#ifndef CUNIMOLSYS_UTILITY_H
#define CUNIMOLSYS_UTILITY_H
#include "resources.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
using std::vector;
using std::string;
using std::pair;


static const double MY_PI = 3.14159265358979323846264338327950;
// For floating point number comparisons
static const double MY_LARGE = 1e9; // Considered to be infinity if more than this.
static const double MY_SMALL = 1e-5; // Considered to be 0 if less than this
/* if a system contains this many atoms, it will be considered as large. In some functions, progress bars will
 * be displayed while processing data */
static const int MY_LARGE_SYSTEM = 200*1000;
/* To avoid memory explosion while accidentally read in a super large file. Limit to 1 billion lines (~10G memory if)
 * each line on average consists of 10 characters. Note: INT_MAX is 2.15 billion */
static const int MY_FILE_LINES_UPPER_BOUND = 1000*1000*1000;

// Resource Files Path
static string DATAFILESPATH=_DATAFILESPATH_;

// String Manipulations

// Split the string according to the separator, returns a vector of string
vector<string> StringSplit(string str,char delimitator=' ');
// Strip away the leading and pending white spaces
string StringStrip(string str);
/* Finds the first occurrence of the token in the string and split the string into a tuple (first,second) without
# the given token. For example StringTok("xyz = 543",'=') returns ("xyz" and "543"). (Both parts are with
# beginning and ending blank spaces striped away. If no such token is found, it returns pair<"","">. */
pair<string,string> StringTok(string str, string token);
string StringToUpper(string str);
string StringToLower(string str);
string StringToCapitalized(string str);
string StringRemoveComment(string str);
bool StringRegexMatch(string str,string pattern,bool caseSensitive=true);
bool StringStartsWith(string str,string pattern,bool caseSensitive=true);
bool StringEndsWith(string str,string pattern,bool caseSensitive=true);

void ProgressBar(double percent,int length = 50);


// This class controls how to deal with stdout and stderr messages.
// It has two child classes: ErrorHandler, OutputHandler
// In a possible GUI version, these messages can be redirected to the Window.
class ErrorAndOutputHandler{
public:
    ErrorAndOutputHandler();
    ~ErrorAndOutputHandler() = default;
    void SetOutput(std::ostream &ostream);
    void TurnOn();
    void TurnOff();
protected:
    bool onOffState_;
    std::ostream *pOutput_;
};

// ErrorHandler,
class ErrorHandler:public ErrorAndOutputHandler{
public:
    ErrorHandler():ErrorAndOutputHandler(){}
    void operator()(string msg,const char* file,int line,bool fatal=true);
};
// OutputHandler,
class OutputHandler:public ErrorAndOutputHandler{
public:
    OutputHandler():ErrorAndOutputHandler(){}
    void operator()(string msg,bool newline=true);
};
// The Singleton Objects
static ErrorHandler error;
// Two macros to quickly output error/warning messages
#define ERROR(msg) error((msg),__FILE__,__LINE__,true)
#define WARNING(msg) error((msg),__FILE__,__LINE__,false)

static OutputHandler output;

// Singleton Object, Initialized only once. User can only access through static members
class PeriodicTable{
public:
    static string AtomicNumberToElement(int atomicNumber);
    static int ElementToAtomicNumber(string element);
    static double AtomicWeight(string element);
    static string PossibleElementWithGivenWeight(double weight);
protected:
    PeriodicTable();
    ~PeriodicTable() = default;
    static PeriodicTable* getInstance_();
    vector<string> elements_;
    std::map<string,double> atomicWeights_;
};

// Mathematical Functions
template <class T=double>
class XYZ_T_{
public:
    XYZ_T_();
    XYZ_T_(T x, T y, T z);
    XYZ_T_(const T xyzs[3]);
    XYZ_T_(const XYZ_T_<T>& otherXYZ) = default;
    T& operator[](int index);
    inline T* ptr() {return xyz;}
    XYZ_T_<T>& operator=(const XYZ_T_<T>& otherXYZ);
    XYZ_T_<T>& operator=(const T xyzs[3]);
    double NormSquared();
    double Norm();
    XYZ_T_<T>& Normalized();
    XYZ_T_<T>& operator+=(XYZ_T_<T> v);
    XYZ_T_<T>& operator-=(XYZ_T_<T> v);
    XYZ_T_<T>& operator*=(double a);
    template<class S>
    friend XYZ_T_<S> operator+(const XYZ_T_<S>& v1, const XYZ_T_<S>& v2);
    template<class S>
    friend XYZ_T_<S> operator-(const XYZ_T_<S>& v1, const XYZ_T_<S>& v2);
    template<class S>
    friend XYZ_T_<S> operator*(const XYZ_T_<S>& v, double a);
    template<class S>
    friend S XYZDot(const XYZ_T_<S> &v1, const XYZ_T_<S> &v2);
    template<class S>
    friend XYZ_T_<S> XYZCross(const XYZ_T_<S> &v1, const XYZ_T_<S> &v2);
    template<class S>
    friend bool operator==(const XYZ_T_<S>& v1, const XYZ_T_<S>& v2);
    template<class S>
    friend std::ostream &operator<<(std::ostream &os, const XYZ_T_<S> &xyz);
private:
    T xyz[3];
};

// Seems like template class must be implemented in .h files. This is quite ugly but seems no other way...
template<class T> XYZ_T_<T>::XYZ_T_():xyz{0, 0, 0}{}
template<class T> XYZ_T_<T>::XYZ_T_(T x, T y, T z):xyz{x, y, z}{}
template<class T> XYZ_T_<T>::XYZ_T_(const T xyzs[3]):xyz{xyzs[0], xyzs[1], xyzs[2]}{}
template<class T> T& XYZ_T_<T>::operator[](int index) {
    if(index<0 || index>2)
        throw std::out_of_range("XYZ[]");
    return xyz[index];
}
template<class T> XYZ_T_<T>& XYZ_T_<T>::operator=(const XYZ_T_<T>& otherXYZ){
    for(int i=0;i<3;i++)
        xyz[i] = otherXYZ.xyz[i];
    return *this;
}
template<class T> XYZ_T_<T>& XYZ_T_<T>::operator=(const T xyzs[3]){
    for(int i=0;i<3;i++)
        xyz[i] = xyzs[i];
    return *this;
}
template<class T> double XYZ_T_<T>::NormSquared(){
    return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2];
}
template<class T> double XYZ_T_<T>::Norm(){
    return sqrt(NormSquared());
}
template<class T> XYZ_T_<T>& XYZ_T_<T>::Normalized(){
    double norm = Norm();
    if(norm < 1e-10)
        throw std::runtime_error("Norm too small: "+ std::to_string(norm));
    (*this)*=1.0/norm;
    return *this;
}
template<class T> XYZ_T_<T>& XYZ_T_<T>::operator+=(XYZ_T_<T> v){
    xyz[0] += v[0];
    xyz[1] += v[1];
    xyz[2] += v[2];
    return *this;
}
template<class T> XYZ_T_<T>& XYZ_T_<T>::operator-=(XYZ_T_<T> v){
    xyz[0] -= v[0];
    xyz[1] -= v[1];
    xyz[2] -= v[2];
    return *this;
}
template<class T> XYZ_T_<T>& XYZ_T_<T>::operator*=(double a){
    xyz[0] *= a;
    xyz[1] *= a;
    xyz[2] *= a;
    return *this;
}
template<class T> XYZ_T_<T> operator+(const XYZ_T_<T>& v1, const XYZ_T_<T>& v2){
    return XYZ_T_<T>(v1.xyz[0] + v2.xyz[0], v1.xyz[1] + v2.xyz[1], v1.xyz[2] + v2.xyz[2]);
}
template<class T> XYZ_T_<T> operator-(const XYZ_T_<T>& v1, const XYZ_T_<T>& v2){
    return XYZ_T_<T>(v1.xyz[0] - v2.xyz[0], v1.xyz[1] - v2.xyz[1], v1.xyz[2] - v2.xyz[2]);
}
template<class T> XYZ_T_<T> operator*(const XYZ_T_<T>& v, double a){
    return XYZ_T_<T>(v.xyz[0] * a, v.xyz[1] * a, v.xyz[2] * a);
}
template<class T> T XYZDot(const XYZ_T_<T> &v1, const XYZ_T_<T> &v2){
    return v1.xyz[0]*v2.xyz[0] + v1.xyz[1]*v2.xyz[1] + v1.xyz[2]*v2.xyz[2];
}
template<class T> XYZ_T_<T> XYZCross(const XYZ_T_<T> &v1, const XYZ_T_<T> &v2){
    double x =  v1.xyz[1]*v2.xyz[2] - v1.xyz[2]*v2.xyz[1];
    double y = -v1.xyz[0]*v2.xyz[2] + v1.xyz[2]*v2.xyz[0];
    double z =  v1.xyz[0]*v2.xyz[1] - v1.xyz[1]*v2.xyz[0];
    return XYZ_T_<T>(x, y, z);
}
template<class T> bool operator==(const XYZ_T_<T>& v1, const XYZ_T_<T>& v2){
    static double tolerance = 1e-5;
    return (v1-v2).Norm() < tolerance;
}
template<class T> std::ostream &operator<<(std::ostream &os, const XYZ_T_<T> &xyz) {
    os <<"["<<xyz.xyz[0]<<", "<<xyz.xyz[1]<<", "<<xyz.xyz[2]<<"]";
    return os;
}
template<class T> void XYZRotate(XYZ_T_<T>* coords, int numAtoms, double degree_clockwise, XYZ_T_<T> axis){
    axis.Normalized();
    double rx = axis[0];
    double ry = axis[1];
    double rz = axis[2];
    // The rotational matrix was copied from a previous code. I forgot its reference...
    double phi = degree_clockwise * MY_PI / 180.0;
    double c = std::cos(phi);
    double s = std::sin(phi);
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
    XYZ_T_<T>* newCoords = new XYZ_T_<T>[numAtoms];
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

using XYZ =  XYZ_T_<double>;
#endif //CUNIMOLSYS_UTILITY_H
