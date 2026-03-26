#include <iostream>
#include <fstream>

#include "methods.hpp"

int main() {
    int steps = 3e4;

    System* systems[] = {
        // new Euler(steps),
        // new EulerSwapped(steps),
        // // new Leapfrog(steps),
        new RK2(steps, 1.0),
        // new RK2(steps, 0.5),
        new RK4(steps),
        // new LinearMultistep(steps, 1)
    };

    std::ofstream f; 
    f.open("index.txt", std::ios::out);

    for (auto &sys : systems) {
        sys->run();
        f << sys->name << std::endl;
    }

    f.close();
}