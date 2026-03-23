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

    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;

    // forces
    double dvx = gamma * m1 * m2 * dx / (dist_squared * dist);
    double dvy = gamma * m1 * m2 * dy / (dist_squared * dist);
    double dvz = gamma * m1 * m2 * dz / (dist_squared * dist);

    // speeds
    v1.x += dvx * dt;
    v1.y += dvy * dt;
    v1.z += dvz * dt;

    v2.x -= dvx * dt;
    v2.y -= dvy * dt;
    v2.z -= dvz * dt;

    // pos
    p1.x += v1.x * dt / m1;
    p1.y += v1.y * dt / m1;
    p1.z += v1.z * dt / m1;

    p2.x += v2.x * dt / m2;
    p2.y += v2.y * dt / m2;
    p2.z += v2.z * dt / m2;

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