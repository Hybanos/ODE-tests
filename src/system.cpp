#include "system.hpp"

void System::run() {

    std::ofstream f;
    f.open(name + ".txt", std::ios::out);

    print();
    save(f);
    while (steps < max_t) {
        t += dt;
        steps += 1;
        step();
        print();
        save(f);
    }

    f.close();
}

void System::print() {
    std::cout << steps << ";"
              << m1 << ";" 
              << p1.x << ";" << p1.y << ";" << p1.z << ";"
              << v1.x << ";" << v1.y << ";" << v1.z << ";"
              << m2 << ";" 
              << p2.x << ";" << p2.y << ";" << p2.z << ";"
              << v2.x << ";" << v2.y << ";" << v2.z 
              << std::endl;
}

void System::save(std::ofstream &f) {
    f << steps << ";"
      << m1 << ";" 
      << p1.x << ";" << p1.y << ";" << p1.z << ";"
      << v1.x << ";" << v1.y << ";" << v1.z << ";"
      << m2 << ";" 
      << p2.x << ";" << p2.y << ";" << p2.z << ";"
      << v2.x << ";" << v2.y << ";" << v2.z 
      << std::endl;
}