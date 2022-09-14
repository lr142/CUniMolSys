#ifndef CUNIMOLSYS_UTILITY_H
#define CUNIMOLSYS_UTILITY_H
#include "resources.h"
#include <iostream>
#include <vector>
using std::vector;
using std::string;
using std::pair;

// Resource Files Path
static string DATAFILESPATH=_DATAFILESPATH_;

// String Manipulation

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


// Singleton Object, Lazy Initialization. User can only access through static members
class PeriodicTable{
public:
    static string AtomicNumberToElement(int atomicNumber);
    static int ElementToAtomicNumber(string element);
    static double AtomicWeight(string element);
protected:
    PeriodicTable();
    ~PeriodicTable() = default;
    static PeriodicTable* getInstance_();
    vector<string> elements_;
    vector<double> atomicWeights_;
};

//class PeriodicTable:
//    def __init__(self):
//        # M is the dummy atom
//        self.elements = "M H He "\
//        "Li Be B C N O F Ne "\
//        "Na Mg Al Si P S Cl Ar "\
//        "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "\
//        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "\
//        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "\
//        "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo";
//        self.elements = self.elements.split()
//        self.atomicweights = {}
//        with open(os.path.join(DATAFILESPATH,"atomicweights.csv"),"r") as awfile:
//            lines = awfile.readlines()
//            for l in lines:
//                words = l.strip().split(",")
//                ele = words[0].upper()
//                weight = float(words[1])
//                self.atomicweights[ele] = weight
//
//    def AtomicNumberToElement(self,i):
//        if i >= len(self.elements):
//            i = 0
//        return self.elements[i]
//    def ElementToAtomicNumber(self,element):
//        for i in range(len(self.elements)):
//            if element == self.elements[i]:
//                return i
//        return -1
//    def AtomicWeight(self,element:str):
//        if element.upper() in self.atomicweights:
//            return self.atomicweights[element.upper()]
//        else:
//            return 0

void f();
int add(int a,int b);
#endif //CUNIMOLSYS_UTILITY_H
