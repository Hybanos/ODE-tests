#include "system.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::time_point;

void System::run() {

    std::ofstream f;
    f.open(name + ".txt", std::ios::out);
    f << "step;m1;p1.x;p1.y;p1.z;v1.x;v1.y;v1.z;m2;p2.x;p2.y;p2.z;v2.x;v2.y;v2.z;K;U;K+1" << std::endl;

    time_point<high_resolution_clock> t1 = high_resolution_clock::now();

    compute_energies();
    save(f);
    while (steps < max_t) {
        t += dt;
        steps += 1;
        step();
        compute_energies();
        save(f);
    }

    f.close();

    time_point<high_resolution_clock> t2 = high_resolution_clock::now();
    std::cout << name << ": ";
    std::cout << (t2 - t1).count() / 1e9 << "s" << std::endl;
}

void System::compute_energies() {
    K = v1.norm() * v1.norm() * m1 / 2 +
        v2.norm() * v2.norm() * m2 / 2;
            
    U = -gamma * (m1 * m2) / (p1 - p2).norm();
}

void System::save(std::ofstream &f) {
    f << steps << ";"
      << m1 << ";" 
      << p1.x << ";" << p1.y << ";" << p1.z << ";"
      << v1.x << ";" << v1.y << ";" << v1.z << ";"
      << m2 << ";" 
      << p2.x << ";" << p2.y << ";" << p2.z << ";"
      << v2.x << ";" << v2.y << ";" << v2.z << ";"
      << K << ";" << U << ";" << K+U
      << std::endl;
}