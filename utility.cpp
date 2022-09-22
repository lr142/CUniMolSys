#include "utility.h"
#include <ctype.h>
#include <fstream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>
#include <regex>
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
bool StringRegexMatch(string str,string pattern,bool caseSensitive){
    regex r;
    if(caseSensitive)
        r = regex(pattern);
    else
        r = regex(pattern,regex::icase);
    smatch sresults;
    if(regex_search(str,sresults,r)){
        return true;
    }
    return false;
}
bool StringStartsWith(string str,string pattern,bool caseSensitive){
    return StringRegexMatch(str,"^"+pattern,caseSensitive);
}
bool StringEndsWith(string str,string pattern,bool caseSensitive){
    return StringRegexMatch(str,pattern+"$",caseSensitive);
}
void ProgressBar(double percent,int length){
    double filled = round(length * percent);
    ostringstream oss;
    oss<<"[";
    int i;
    for(i=0;i<filled;i++)
        oss<<"#";
    for(i=filled;i<length;i++){
        oss<<" ";
    }
    oss<<"] "<<fixed<<setprecision(2)<<setw(4)<<percent*100<<"%";
    if(percent > 1.0-MY_SMALL)
        oss<<"\n";
    cout<<"\r"+oss.str()<<ends<<flush;
}

/* These two functions will be convenient when reading files if you want to jump to a desired
 * location. The 1st version scans lines between line lo and line hi for pattern. If found, it
 * writes the correct line no to lineno and returns true. If not found, it returns false.
 * In the 2nd version, an fstream is read at most max_lines lines. If found, the desired line
 * is written in line and returns true, otherwise returns false. */
bool JumpToLine(vector<string> &lines, string pattern, int &lineno, int lo, int hi){
    for(lineno=lo;lineno<hi;lineno++){
        if(StringRegexMatch(lines[lineno],pattern))
            return true;
    }
    return false;
}
bool JumpToLine(ifstream &ifs, string pattern, string &line, int max_lines){
    int counter = 0;
    while( (not ifs) and counter<max_lines ){
        counter++;
        getline(ifs,line);
        if(StringRegexMatch(line,pattern))
            return true;
    }
    return false;
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

// Multi-threading related functions
pair<int,int> TaskDistribution(int iThread,int nThreads,int nTasks){
    int tasks_per_thread = nTasks%nThreads==0? nTasks/nThreads : nTasks/nThreads+1 ;
    int start = tasks_per_thread*iThread;
    int end = min(tasks_per_thread*(iThread+1),nTasks);
    return make_pair(start,end);
}