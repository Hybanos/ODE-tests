#include <iostream>
#include <fstream>

#include "system.hpp"

template <typename T>
void haha(config c, bool save) {
    T s(c, save);
    s.run();
}

void run_all(config c, bool save) {
    haha<System<Euler>>(c, save);
    haha<System<EulerSymplectic>>(c, save);
    haha<System<LeapFrog>>(c, save);
    haha<System<RK2>>(c, save);
    haha<System<RK4>>(c, save);
    haha<System<RK45>>(c, save);
    haha<System<DOP853>>(c, save);
}

int main(int argc, char *argv[]) {
    int bodies;
    double target_t, seed;
    bool save = false;

    if (argc < 4) exit(1);
    bodies = std::stoi(argv[1]);
    target_t = std::stod(argv[2]);
    seed = std::stod(argv[3]);
    if (argc == 5 && argv[4][0] == 't') save = true;

    config c(bodies, target_t, seed);
    run_all(c, save);

    return 0;
}
