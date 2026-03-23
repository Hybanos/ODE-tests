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

    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;

    // forces
    double dvx1 = -gamma * m2 * dx / dist_squared * dist;
    double dvy1 = -gamma * m2 * dy / dist_squared * dist;
    double dvz1 = -gamma * m2 * dz / dist_squared * dist;

    double dvx2 = gamma * m1 * dx / dist_squared * dist;
    double dvy2 = gamma * m1 * dy / dist_squared * dist;
    double dvz2 = gamma * m1 * dz / dist_squared * dist;

    // speeds
    v1.x += dvx1 * dt;
    v1.y += dvy1 * dt;
    v1.z += dvz1 * dt;

    v2.x += dvx2 * dt;
    v2.y += dvy2 * dt;
    v2.z += dvz2 * dt;

    // pos
    p1.x += v1.x;
    p1.y += v1.y;
    p1.z += v1.z;

    p2.x += v2.x;
    p2.y += v2.y;
    p2.z += v2.z;

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