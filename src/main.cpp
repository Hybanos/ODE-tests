#include <iostream>

#include "methods.hpp"

int main() {
    int steps = 30000;

    Euler s1(steps);
    s1.run();

    // RK21
    RK2 s2(steps, 1.0);
    s2.run();

    // RK22
    RK2 s3(steps, 0.5);
    s3.run();
}