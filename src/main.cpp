#include <iostream>
#include <fstream>

#include "methods.hpp"

int main() {
    int steps = 3e4;
    int bodies = 4;
    int seed = 40;

    System* systems[] = {
        // new Euler(steps, bodies, seed),
        new EulerSwapped(steps, bodies, seed),
        // new Leapfrog(steps, bodies, seed),
        // new RK2(steps, bodies, seed, 1.0),
        // new RK2(steps, bodies, seed, 0.5),
        // new RK4(steps, bodies, seed),
        new LinearMultistep(steps, 1, bodies, seed)
    };

    std::ofstream f; 
    f.open("index.txt", std::ios::out);

    for (auto &sys : systems) {
        sys->run();
        f << sys->name << std::endl;
    }

    f.close();
}