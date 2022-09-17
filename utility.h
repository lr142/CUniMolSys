#ifndef CUNIMOLSYS_UTILITY_H
#define CUNIMOLSYS_UTILITY_H
#include "resources.h"
#include <iostream>
#include <vector>
#include <map>
using std::vector;
using std::string;
using std::pair;

// Resource Files Path
static string DATAFILESPATH=_DATAFILESPATH_;

// String Manipulations

// Split the string according to the separator, returns a vector of string
vector<string> StringSplit(string str,char delimitator=' ');
// Strip away the leading and pending white spaces
string StringStrip(string str);
// Test whether the string starts with the pattern (case sensitive). If the pattern is empty, it will return true
bool StringStartsWith(string str,string pattern);
/* Finds the first occurrence of the token in the string and split the string into a tuple (first,second) without
# the given token. For example StringTok("xyz = 543",'=') returns ("xyz" and "543"). (Both parts are with
# beginning and ending blank spaces striped away. If no such token is found, it returns pair<"","">. */
pair<string,string> StringTok(string str, string token);
string StringToUpper(string str);
string StringToLower(string str);
string StringToCapitalized(string str);
string StringRemoveComment(string str);

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
class XYZ{
public:
    XYZ();
    XYZ(double x,double y,double z);
    XYZ(const double xyzs[3]);
    XYZ(const XYZ& otherXYZ) = default;
    double& operator[](int index);
    inline double* ptr() {return xyz;}
    XYZ& operator=(const XYZ& otherXYZ);
    XYZ& operator=(const double xyzs[3]);
    double NormSquared();
    double Norm();
    XYZ& Normalized();
    XYZ& operator+=(XYZ v);
    XYZ& operator-=(XYZ v);
    XYZ& operator*=(double a);
    friend XYZ operator+(const XYZ& v1, const XYZ& v2);
    friend XYZ operator-(const XYZ& v1, const XYZ& v2);
    friend XYZ operator*(const XYZ& v, double a);
    friend double XYZDot(const XYZ &v1, const XYZ &v2);
    friend XYZ XYZCross(const XYZ &v1, const XYZ &v2);
    friend bool operator==(const XYZ& v1, const XYZ& v2);
    friend std::ostream &operator<<(std::ostream &os, const XYZ &xyz);
private:
    double xyz[3];
};
void XYZRotate(XYZ* coords,int numAtoms,double degree_clockwise,XYZ axis);
#endif //CUNIMOLSYS_UTILITY_H
