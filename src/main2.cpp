#include <iostream>
#include <fstream>

#include "system2.hpp"

int main() {
    config c(2, 30);

    System2<RK4> s(c);
    s.run();
}