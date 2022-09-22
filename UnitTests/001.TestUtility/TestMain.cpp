#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include <vector>
#include <fstream>
using namespace std;

TEST(XYZ,template){
    float v[] = {1.1,2.,3.5};
    auto x = XYZ_T_<float>(v);
    cout<<x<<" "<<sizeof(x)<<endl;

    double w[] = {1.1,2.,3.5};
    auto y = XYZ_T_<double>(w);
    cout<<y<<" "<<sizeof(y)<<endl;
}

TEST(StringFunctions,split){
    vector<string> result = {"1","2","3"};
    EXPECT_EQ(StringSplit("1 2  3 "),result);
    EXPECT_EQ(StringSplit("  1 2  3"),result);
}
TEST(StringFunctions,strip){
    EXPECT_EQ(StringStrip("    "),"");
    EXPECT_EQ(StringStrip(" "),"");
    EXPECT_EQ(StringStrip("\t\n"),"");
    EXPECT_EQ(StringStrip("   A"),"A");
    string result = "works fine";
    EXPECT_EQ(StringStrip("  works fine   "),result);
    EXPECT_EQ(StringStrip("\r\nworks fine\r"),result);
}
TEST(StringFunctions,startswith){
    EXPECT_EQ(StringStartsWith("  ABC","  A"),true);
    EXPECT_EQ(StringStartsWith("  ABC",""),true);
    EXPECT_EQ(StringStartsWith("xyz","xyzz"),false);
    EXPECT_EQ(StringStartsWith("x\tyz","x\ty"),true);
}
TEST(StringFunctions,tok){
    pair<string,string> result = make_pair("A","BC");
    EXPECT_EQ(StringTok("  ABC","  A"),result);
    result = make_pair("AC","BC");
    EXPECT_EQ(StringTok("AC BC","  AC    "),result);
    result = make_pair("","");
    EXPECT_EQ(StringTok("ACBC","  AB    "),result);
}
TEST(StringFunctions,other){
    EXPECT_EQ(StringToUpper("  abc"),"  ABC");
    EXPECT_EQ(StringToUpper("1 21AaBBcDea"),"1 21AABBCDEA");
    EXPECT_EQ(StringRemoveComment("1 21AaBBcD#ea"),"1 21AaBBcD");
    EXPECT_EQ(StringRemoveComment("# da"),"");
    EXPECT_EQ(StringRemoveComment(" # da"),"");
    EXPECT_EQ(StringRemoveComment("ab ### da"),"ab");
}
TEST(ErrorHandler,tok){
    WARNING("Show some error message 1");
    error.TurnOff();
    WARNING("Show some error message 2");
    error.TurnOn();
    WARNING("Show some error message 3");
    error.TurnOff();
    ofstream ofs("dump.txt");
    error.SetOutput(ofs);
    WARNING("Show some error message 4");
    error.SetOutput(cout);
    error.TurnOn();
    WARNING("Show some error message 5");

    output("Show some output message 1");
    output.TurnOff();
    output("Show some output message 2");
    output.TurnOn();
    output("Show some output message 3");
    output.TurnOff();
    output.SetOutput(ofs);
    output("Show some output message 4");
    output.SetOutput(cout);
    output.TurnOn();
    output("Show some output message 5");
}

TEST(TestPeriodicTable,case1){
    EXPECT_EQ(PeriodicTable::AtomicNumberToElement(1),"H");
    EXPECT_EQ(PeriodicTable::AtomicNumberToElement(10),"Ne");
    EXPECT_EQ(PeriodicTable::ElementToAtomicNumber("O"),8);
    EXPECT_EQ(PeriodicTable::ElementToAtomicNumber("LR"),103);
    EXPECT_NEAR(PeriodicTable::AtomicWeight("LR"),266,1e-4);
    EXPECT_NEAR(PeriodicTable::AtomicWeight("pT"),195.08,1e-4);
    EXPECT_NEAR(PeriodicTable::AtomicWeight("M"),1e-30,1e-4);
    EXPECT_EQ(PeriodicTable::PossibleElementWithGivenWeight(23.14),"Na");
    EXPECT_EQ(PeriodicTable::PossibleElementWithGivenWeight(12.3),"C");
    EXPECT_EQ(PeriodicTable::PossibleElementWithGivenWeight(0),"M");
}

TEST(XYZ,case1){
    XYZ v1 = {1.0,2,3};
    XYZ v2(4,5,6.0);
    XYZ_DTYPE coords[] = {5,7,9};
    XYZ_DTYPE *pD = coords;
    XYZ v3 = pD;
    EXPECT_EQ(v1[0],1.0);
    EXPECT_EQ(v1[1],2.0);
    EXPECT_EQ(v1[2],3.0);
    EXPECT_TRUE(v1+v2==v3);
    EXPECT_EQ(v1-v2,XYZ(-3,-3,-3));
    EXPECT_NEAR(XYZDot(v1, v2), 32.0, 1e-6);
    EXPECT_TRUE(v3*1.1==XYZ({5.5,7.7,9.9}));

    XYZ v4 = {1.5, 2.9, -0.8};
    EXPECT_TRUE(v4.Normalized()==XYZ(0.446223, 0.862698, -0.237986));
    //cout<<v4.Normalized()<<endl;
    v4 = {1.5, 2.9, -0.8};
    v4-={0.5,0.4,1.2};
    EXPECT_TRUE(v4==XYZ(1.0,2.5,-2.0));
    v4+={0.5,0.4,1.2};
    EXPECT_TRUE(v4==XYZ(1.5, 2.9, -0.8));

    XYZ v5 = {0,0,0.0};
    EXPECT_ANY_THROW(v5.Normalized());

    XYZ v6 = {1.5,6.2,-9.4};
    XYZ v7 = {2,8.5,-1};
    XYZ v8 = {73.7,-17.3,0.35};
    EXPECT_TRUE(XYZCross(v6, v7) == v8);
}
TEST(XYZ,rotate){
    XYZ coords[4];
    coords[0] = {1,0,0};
    coords[1] = {0,1,0};
    coords[2] = {0,0,1};
    coords[3] = {1,1,1};
    XYZ direction = {0,0,1};
    XYZRotate(coords,4,120,{100,100,100});
    // x-->z, y-->x, z-->y,  [111] is unchanged.
    EXPECT_TRUE(coords[0]==XYZ(0,0,1));
    EXPECT_TRUE(coords[1]==XYZ(1,0,0));
    EXPECT_TRUE(coords[2]==XYZ(0,1,0));
    EXPECT_TRUE(coords[3]==XYZ(1,1,1));
}