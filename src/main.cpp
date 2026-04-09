#include <iostream>
#include <fstream>

#include "system.hpp"

config c(4, 40);
std::ofstream f; 

template <typename T>
void haha() {
    T s(c);
    s.run();
    f << s.get_name() << std::endl;
}

int main() {
    f.open("index.txt", std::ios::out);

    // haha<System2<Euler>>();
    // haha<System2<EulerSymplectic>>();
    // haha<System2<LeapFrog>>();
    // haha<System2<RK2>>();
    // haha<System2<RK4>>();
    haha<System<RK45>>();
    haha<System<DOP853>>();

    f.close();
}
