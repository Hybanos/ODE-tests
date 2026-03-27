#include "system.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::time_point;


System::System(std::string _name, int _max_t, int bodies = 2, int seed = 0) : name{_name}, max_t{_max_t} {
    m.resize(bodies);
    x.resize(bodies);
    v.resize(bodies);

    // generate random masses/positions/speeds
    std::srand(seed);
    for (int i = 0; i < bodies; i++) {
        m[i] = ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2;
        x[i] = vec3{
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2,
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2,
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2,
        };
        v[i] = vec3{
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2,
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2,
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 2,
        };
    }

    // offset those values so the barycenter is the origin of the coordinates
    vec3 barycenter_speed = vec3{0, 0, 0};
    vec3 barycenter_pos = vec3{0, 0, 0};
    double m_sum = 0;
    for (int i = 0; i < bodies; i++) {
        barycenter_speed = barycenter_speed + v[i] * m[i];
        barycenter_pos = barycenter_pos + x[i] * m[i];
        m_sum += m[i];
    }
    barycenter_speed = barycenter_speed / m_sum;
    barycenter_pos = barycenter_pos / m_sum;

    for (int i = 0; i < bodies; i++) {
        v[i] = v[i] - barycenter_speed;
        x[i] = x[i] - barycenter_pos;
    }
};

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
        if (!(steps % 1)) save(f);
    }

    f.close();

    time_point<high_resolution_clock> t2 = high_resolution_clock::now();
    std::cout << name << ": ";
    std::cout << (t2 - t1).count() / 1e9 << "s" << std::endl;
}

vec3 System::compute_acceleration(vec3 x1, vec3 x2) {
    double dist_squared = x1.dist_squared(x2);
    double dist = std::sqrt(dist_squared);

    vec3 r  = x2 - x1;
    vec3 dv = - r * gamma / (dist_squared * dist);

    return dv;
}

void System::compute_energies() {
    double m1 = m[0];
    double m2 = m[1];
    vec3 p1 = x[0];
    vec3 p2 = x[1];
    vec3 v1 = v[0];
    vec3 v2 = v[1];

    K = v1.norm() * v1.norm() * m1 / 2 +
        v2.norm() * v2.norm() * m2 / 2;
            
    U = -gamma * (m1 * m2) / (p1 - p2).norm();
}

void System::save(std::ofstream &f) {
    double m1 = m[0];
    double m2 = m[1];
    vec3 p1 = x[0];
    vec3 p2 = x[1];
    vec3 v1 = v[0];
    vec3 v2 = v[1];

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