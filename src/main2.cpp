#include <iostream>
#include <fstream>

#include "system2.hpp"

config c(2, 30);
std::ofstream f; 

template <typename T>
void haha() {
    T s(c);
    s.run();
    f << s.get_name() << std::endl;
}

int main() {
    f.open("index.txt", std::ios::out);

    haha<System2<Euler>>();
    haha<System2<RK4>>();

    f.close();
}
