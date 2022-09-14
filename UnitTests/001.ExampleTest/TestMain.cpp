#include <iostream>
#include "gtest/gtest.h"
#include "utility.h"

TEST(Example,t1){
    EXPECT_EQ(add(1,2),3);
}
TEST(Example,t2){
    EXPECT_EQ(add(3,4),7);
}
TEST(Example,t3){
    EXPECT_EQ(add(11,12),22);
}
