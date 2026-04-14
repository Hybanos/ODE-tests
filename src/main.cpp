#include <iostream>
#include <fstream>

#include "system.hpp"

config c(200, 10, 30);
std::ofstream f; 

template <typename T>
void haha() {
    T s(c);
    s.run();
    f << s.get_name() << std::endl;
}

int main() {
    f.open("index.txt", std::ios::out);

    // haha<System<Euler>>();
    // haha<System<EulerSymplectic>>();
    // // haha<System<LeapFrog>>();
    haha<System<RK2>>();
    // haha<System<RK4>>();
    // haha<System<RK45>>();
    // haha<System<DOP853>>();

    f.close();
}
