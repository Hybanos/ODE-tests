#include "system.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::time_point;


System::System(std::string _name, double _target_t, int bodies = 2, int seed = 0) : name{_name}, target_t{_target_t} {
    m.resize(bodies);
    x.resize(bodies);
    v.resize(bodies);
    a.resize(bodies);

    // generate random masses/positions/speeds
    std::srand(seed);
    for (int i = 0; i < bodies; i++) {
        m[i] = ((double) (std::rand()) / RAND_MAX) * 2 + 3;
        x[i] = vec3{
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 4,
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 4,
            ((double) (std::rand()) / RAND_MAX * 2 - 1) * 4,
        };
        v[i] = vec3{
            ((double) (std::rand()) / RAND_MAX * 2 - 1),
            ((double) (std::rand()) / RAND_MAX * 2 - 1),
            ((double) (std::rand()) / RAND_MAX * 2 - 1),
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
    f << "nbody;step;t;";
    for (int i = 0; i < m.size(); i++) {
        f << "m" << i << ";";
        f << "x" << i << ".x;";
        f << "x" << i << ".y;";
        f << "x" << i << ".z;";
        f << "v" << i << ".x;";
        f << "v" << i << ".y;";
        f << "v" << i << ".z;";
    }
    f << "K;U;K+U" << std::endl;;

    time_point<high_resolution_clock> t1 = high_resolution_clock::now();

    compute_energies();
    save(f);
    while (t < target_t) {
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

void System::compute_accelerations(std::vector<vec3> &a, std::vector<vec3> &x) {
    for (int i = 0; i < m.size(); i++) {
        vec3 a_i = vec3{0, 0, 0};

        for (int j = 0; j < i; j++) {
            vec3 r = x[j] - x[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        for (int j = i+1; j < m.size(); j++) {
            vec3 r = x[j] - x[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }
        a[i] = a_i;
    }
}

void System::compute_energies() {
    double m1 = m[0];
    double m2 = m[1];
    vec3 p1 = x[0];
    vec3 p2 = x[1];
    vec3 v1 = v[0];
    vec3 v2 = v[1];

    K = 0;
    U = 0;
    for (int i = 0; i < m.size(); i++) {
        K += v[i].norm() * v[i].norm() * m[i] / 2;
        for (int j = i+1; j < m.size(); j++) {
            U += -gamma * m[i] * m[j] / (x[j] - x[i]).norm();
        }
    }

    // U = -gamma * (m1 * m2) / (p1 - p2).norm();
}

void System::save(std::ofstream &f) {
    double m1 = m[0];
    double m2 = m[1];
    vec3 p1 = x[0];
    vec3 p2 = x[1];
    vec3 v1 = v[0];
    vec3 v2 = v[1];

    f << m.size() << ";" << steps << ";" << t << ";";
    for (int i = 0; i < m.size(); i++) {
        f << m[i] << ";" 
          << x[i].x << ";" << x[i].y << ";" << x[i].z << ";"
          << v[i].x << ";" << v[i].y << ";" << v[i].z << ";";
    } 

    f << K << ";" << U << ";" << K+U
      << std::endl;
}