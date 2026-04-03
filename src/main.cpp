#include <iostream>
#include <fstream>

#include "methods.hpp"

int main() {
    double target_t = 30;
    int bodies = 2;
    int seed = 30;

    System* systems[] = {
        // new Exact(target_t, bodies, seed),
        // new Euler(target_t, bodies, seed),
        // new EulerSwapped(target_t, bodies, seed),
        // new Leapfrog(target_t, bodies, seed),
        // new RK2(target_t, bodies, seed, 1.0),
        // new RK2(target_t, bodies, seed, 0.5),
        new RK4(target_t, bodies, seed),
        // new RK45(target_t, bodies, seed),
        new DOP853(target_t, bodies, seed),
        // new LinearMultistep(target_t, 1, bodies, seed)
    };

    std::ofstream f; 
    f.open("index.txt", std::ios::out);

    for (auto &sys : systems) {
        sys->run();
        f << sys->name << std::endl;
    }

    f.close();
}