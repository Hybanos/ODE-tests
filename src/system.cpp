#include "system.hpp"

void System::run() {
    print();
    while (t < max_t) {
        t += dt;
        step();
    }
}

void System::step() {
    double dist_squared = p1.dist_squared(p2);
    double dist = std::sqrt(dist_squared);

    vec3 dx = p2 - p1;

    // forces
    vec3 dv = dx * gamma * m1 * m2 / (dist_squared * dist);

    // speeds
    v1 = v1 + dv * dt;
    v2 = v2 - dv * dt;

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