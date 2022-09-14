#include "utility.h"
#include <ctype.h>
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


// Singleton Object, Lazy Initialization. User can only access through static members
//class PeriodicTable{
//public:
//    static std::string AtomicNumberToElement(int atomicNumber);
//    static int ElementToAtomicNumber(std::string element);
//    static double AtomicWeight(std::string element)
//protected:
//    PeriodicTable();
//    ~PeriodicTable() = default;
//    static PeriodicTable* getInstance_();
//    std::vector<string> elements_;
//    std::vector<double> atomicWeights_;
//};
PeriodicTable::PeriodicTable() {
    elements_ = StringSplit("M H He "\
        "Li Be B C N O F Ne "\
        "Na Mg Al Si P S Cl Ar "\
        "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "\
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "\
        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "\
        "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo");
}