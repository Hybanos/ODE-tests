#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>

#include <mdspan/mdspan.hpp>

#include "vec3.hpp"
#include "odes.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::time_point;

struct config {
    int bodies;
    double target_t;

    std::vector<vec3> x;
    std::vector<vec3> v;
    std::vector<double> m;

    config(int _bodies, int _seed);
};

template <typename Solver>
class System2 {
    private:
        // initial configuration
        config c;

        Solver *solver = nullptr;

        int bodies;
        double gamma = 1;
        // double target_t = 0.001;
        double target_t = 30;

        std::vector<double> _data;
        array data;

        // kinetic and potential energies
        double K, U;

        void compute_ac(double t, array &Y, array &ret);
        void compute_energies();
        void save(std::ofstream &f);
    public:
        void run();
        System2(config &_c);
        ~System2() {delete solver;}
};

template <typename Solver>
System2<Solver>::System2(config &_c) : c{_c} {
    bodies = c.bodies;
    for (int i = 0; i < bodies; i++) {
        _data.push_back(c.x[i].x);
        _data.push_back(c.x[i].y);
        _data.push_back(c.x[i].z);
    }

    for (int i = 0; i < bodies; i++) {
        _data.push_back(c.v[i].x);
        _data.push_back(c.v[i].y);
        _data.push_back(c.v[i].z);
    }

    for (int i = 0; i < bodies; i++) {
        _data.push_back(c.m[i]);
    }

    data = array(_data.data(), bodies * 7);
    solver = new Solver([&](double t, array& Y, array& ret){this->compute_ac(t, Y, ret);}, data);
}

template <typename Solver>
void System2<Solver>::compute_ac(double t, array &Y, array &ret) {
    vecarray x = vecarray((vec3 *) Y.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) Y.data_handle() + bodies, bodies);
    array m = array(Y.data_handle() + bodies * 6, bodies);
    vecarray ret_v = vecarray((vec3 *) ret.data_handle(), bodies);
    vecarray ret_a = vecarray((vec3 *) ret.data_handle() + bodies, bodies);

    for (int i = 0; i < bodies; i++) {
        vec3 a_i = vec3{0, 0, 0};

        for (int j = 0; j < i; j++) {
            vec3 r = x[j] - x[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        for (int j = i+1; j < bodies; j++) {
            vec3 r = x[j] - x[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        ret_a[i] = a_i / m[i];
        ret_v[i] = v[i];
    } 
}

template <typename Solver>
void System2<Solver>::run() {
    std::ofstream f;
    f.open(solver->name + ".txt", std::ios::out);
    f << "nbody;step;t;";
    for (int i = 0; i < bodies; i++) {
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
    while (solver->t < target_t) {
        solver->t += solver->dt;
        solver->step();
        compute_energies();
        save(f);
    }

    f.close();
    std::cout << "haha" << std::endl;

    time_point<high_resolution_clock> t2 = high_resolution_clock::now();
    std::cout << solver->name << ": ";
    std::cout << (t2 - t1).count() / 1e9 << "s" << std::endl;
}

template <typename Solver>
void System2<Solver>::compute_energies() {
    size_t size = bodies * 3;
    vecarray x = vecarray((vec3 *) data.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) data.data_handle() + bodies, bodies);
    array m = array(data.data_handle() + bodies * 6, bodies);

    K = 0;
    U = 0;
    for (int i = 0; i < bodies; i++) {
        K += v[i].norm() * v[i].norm() * m[i] / 2;
        for (int j = i+1; j < m.size(); j++) {
            U += -gamma * m[i] * m[j] / (x[j] - x[i]).norm();
        }
    }
}

template <typename Solver>
void System2<Solver>::save(std::ofstream &f) {
    vecarray x = vecarray((vec3 *) data.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) data.data_handle() + bodies, bodies);
    array m = array(data.data_handle() + bodies * 6, bodies);

    f << bodies << ";" << solver->steps << ";" << solver->t << ";";
    for (int i = 0; i < bodies; i++) {
        f << m[i] << ";" 
          << x[i].x << ";" << x[i].y << ";" << x[i].z << ";"
          << v[i].x << ";" << v[i].y << ";" << v[i].z << ";";
    }

    f << K << ";" << U << ";" << K+U
      << std::endl;
}