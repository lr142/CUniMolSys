#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"
#include <vector>
using namespace std;

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
