#include "system.hpp"

void System::run() {
    print();
    while (t < max_t) {
        t += dt;
        step();
    }
}

void System::step() {

    // speeds
    v1 = v1 + euler(f, dt);
    v2 = v2 - euler(f, dt); 

    // pos
    p1 = p1 + v1 * dt / m1;
    p2 = p2 + v2 * dt / m2;

    print();
}

void System::print() {
    std::cout << t  << ";"
              << m1 << ";" 
              << p1.x << ";" << p1.y << ";" << p1.z << ";"
              << v1.x << ";" << v1.y << ";" << v1.z << ";"
              << m2 << ";" 
              << p2.x << ";" << p2.y << ";" << p2.z << ";"
              << v2.x << ";" << v2.y << ";" << v2.z 
              << std::endl;
}