#include <iostream>
#include <fstream>

#include "system.hpp"

template <typename T>
void haha(config c) {
    T s(c, false);
    s.run();
}

void run_all(config c) {
    haha<System<Euler>>(c);
    haha<System<EulerSymplectic>>(c);
    haha<System<LeapFrog>>(c);
    haha<System<RK2>>(c);
    haha<System<RK4>>(c);
    haha<System<RK45>>(c);
    haha<System<DOP853>>(c);
}

int main(int argc, char *argv[]) {
    int bodies;
    double target_t, seed;
    bool save = false;

    if (argc < 4) exit(1);
    bodies = std::stoi(argv[1]);
    target_t = std::stod(argv[2]);
    seed = std::stod(argv[3]);
    if (argc == 5 && argv[4] == "true") save = true;

    config c(bodies, target_t, seed);
    run_all(c);

    return 0;
}
