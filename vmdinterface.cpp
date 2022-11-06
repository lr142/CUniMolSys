#include "vmdinterface.h"
#include <sstream>
using namespace std;

VMDInterface::VMDInterface(string molecularFileName, std::string scriptFileName, bool autobonds) {
    this->molecularFileName = molecularFileName;
    this->scriptFileName = scriptFileName;
    this->autobonds = autobonds;
    this->textsize = 1.5;
    this->textthickness = 1.5;
    this->commands = vector<string>();
    Step1_Representations();
    Step2_LoadMolecule();
    Step3_Representations();
}

void VMDInterface::Step1_Representations() {
    // Adjust Colors for some elements
    AddCommand("color	Element	C	gray");
    AddCommand("color	Element	Al	magenta");
    AddCommand("color	Element	Si	orange");
    AddCommand("color	Element	Mg	green");
    AddCommand("color	Element	Na	blue");
    AddCommand("color	Element	K	purple");
    AddCommand("color	Element	Cl	green");
    AddCommand("#display cuedensity 0.0");  // Turn off depth cue if you like by uncommenting it
    AddCommand("display projection Orthographic");
    AddCommand("menu main on");
    AddCommand("mol modstyle 0 0 Licorice 0.300000 12.000000 12.000000");
    AddCommand("mol color Name");
    AddCommand("mol representation Licorice 0.300000 12.000000 12.000000");
    AddCommand("mol selection all");
    AddCommand("mol material Opaque");
    AddCommand("mol modrep 0 0");
    AddCommand("color Display Background white");
}

void VMDInterface::Step2_LoadMolecule() {
    // Autobonds On/Off 可选择是否自动加载成键。关闭的话，可用于检测体系中成键是否正常
    ostringstream oss;
    oss<<"mol new "<<molecularFileName<<" autobonds "<< (autobonds?"on":"off");
    AddCommand(oss.str());
    if(StringEndsWith(molecularFileName,"mol2",false)){
        AddCommand("set sel [atomselect top all]");
        AddCommand("$sel set element [$sel get type]");
    }
}

void VMDInterface::Step3_Representations() {
    AddCommand("mol modcolor 0 0 Element");
    // 以下根据原子数判断。如果 > 1000, 用Lines，否则用 Stick-and-Balls
    AddCommand("lassign [molinfo top get numatoms] NAtoms");
    AddCommand("if {$NAtoms < 1000} {;mol modstyle top 0 Licorice;} else {mol modstyle top 0 Lines}");
    for(string keyword:vector<string>{"Atoms","Bonds", "Angles", "Dihedrals"})
        AddCommand("color Labels "+keyword+" black");
    AddCommand("label textthickness 3.6");
    AddCommand("label textsize 1.0");
    //label add Atoms 0/{i}
    //label textoffset Atoms 0 { -0.007692 -0.060440 } # 设置偏移
    //label textformat Atoms 0 { custem  }  # 设置Label内容。
    //label textformat Atoms 0 { %i  } # %加字符是一些特殊标记，包括：%a name, %q charge, %i 0-index, %1i 1-index, %e element
}

void VMDInterface::AddLabels(vector<int> indexes, vector<string> labels, double *offset) {
    string off;
    if(offset == nullptr)
        off = "{0.03 0.03 0.03}";
    else{
        ostringstream oss;
        oss<<"{"<<offset[0]<<" "<<offset[1]<<" "<<offset[2]<<"}";
        off = oss.str();
    }
    for(int i=0;i<indexes.size();i++){
        string ind = to_string(indexes[i]);
        AddCommand("label add Atoms 0/"+ ind);
        AddCommand("label textoffset Atoms "+ ind +" "+off);
        AddCommand("label textformat Atoms "+ ind + " {"+labels[i]+"}");
    }
}

string __put_6_numbers_in_a_line(double a,double b,double c,double d,double e,double f){
    ostringstream oss;
    oss<<"{"<< a <<" "<< b <<" " << c <<" } { " << d <<" "<< e <<" " << f <<" }";
    return oss.str();
}
void VMDInterface::AddPeriodicBoundaryBox(Boundary boundary, std::string color, int width, std::string style) {
    // Boundary in the [[xlo,xhi],[ylo,yhi],[zlo,zhi]] format
    AddCommand("graphics 0 color "+color);
    XYZ_DTYPE lohi[3][2];
    boundary.GetLoHi(lohi);
    double xlo = lohi[0][0];double xhi = lohi[0][1];
    double ylo = lohi[1][0];double yhi = lohi[1][1];
    double zlo = lohi[2][0];double zhi = lohi[2][1];

    string lines[] = {
            __put_6_numbers_in_a_line(xlo, ylo, zlo, xlo, yhi, zlo),
            __put_6_numbers_in_a_line(xhi, ylo, zlo, xhi, yhi, zlo),
            __put_6_numbers_in_a_line(xlo, ylo, zlo, xhi, ylo, zlo),
            __put_6_numbers_in_a_line(xlo, yhi, zlo, xhi, yhi, zlo), // Bottom plane
            __put_6_numbers_in_a_line(xlo, ylo, zhi, xlo, yhi, zhi),
            __put_6_numbers_in_a_line(xhi, ylo, zhi, xhi, yhi, zhi),
            __put_6_numbers_in_a_line(xlo, ylo, zhi, xhi, ylo, zhi),
            __put_6_numbers_in_a_line(xlo, yhi, zhi, xhi, yhi, zhi), // Top plane
            __put_6_numbers_in_a_line(xlo, ylo, zlo, xlo, ylo, zhi),
            __put_6_numbers_in_a_line(xlo, yhi, zlo, xlo, yhi, zhi),
            __put_6_numbers_in_a_line(xhi, ylo, zlo, xhi, ylo, zhi),
            __put_6_numbers_in_a_line(xhi, yhi, zlo, xhi, yhi, zhi) // Four pillars
    };
    for(int i=0;i<12;i++)
        AddCommand("graphics 0 line "+lines[i]+" width "+ to_string(width)+" style " + style);
}

void VMDInterface::AddCommand(string cmd) {
    this->commands.push_back(cmd);

}

void VMDInterface::Run(){
    ofstream ofs(scriptFileName);
    for(auto &cmd:commands){
        ofs<<cmd<<endl;
    }
}
