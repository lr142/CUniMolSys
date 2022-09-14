#ifndef CUNIMOLSYS_UTILITY_H
#define CUNIMOLSYS_UTILITY_H
#include <iostream>
class PeriodicTable{
public:
    PeriodicTable();
    ~PeriodicTable();
    std::string AtomicNumberToElement(int atomicNumber);
    int ElementToAtomicNumber(std::string element);
    double f();
};

void f();
int add(int a,int b);
#endif //CUNIMOLSYS_UTILITY_H
