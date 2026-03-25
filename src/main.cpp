#include <iostream>

#include "methods.hpp"

int main() {
    int steps = 3000;

    Euler s1(steps);
    s1.run();

    RK21 s2(steps);
    s2.run();
}