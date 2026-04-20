#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>
#include <iomanip>

#include <mdspan/mdspan.hpp>

#include "vec3.hpp"
#include "odes.hpp"
#include "defs.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::time_point;

struct config {
    int bodies;
    int seed;
    fpoint_t target_t;

    std::vector<vec3> x;
    std::vector<vec3> v;
    std::vector<fpoint_t> m;

    config(int _bodies, fpoint_t _target_t, int _seed);
};

template <typename Integrator>
class System {
    private:
        // initial configuration
        config c;

        Integrator *integrator = nullptr;

        int bodies;
        fpoint_t gamma = 1;
        fpoint_t eps_r = 1e-6;
        // fpoint_t eps_r = 0;

        std::vector<fpoint_t> _data;
        array data;
        array x;
        array v;
        array m;

        // kinetic and potential energies
        fpoint_t K, U;
        // angular momentum
        vec3 A;

        bool save_results;

        void compute_ac(fpoint_t t, array &Y, array &ret);
        void compute_ac_order_2(fpoint_t t, array &x, array &ret);
        void compute_energies();
        void compute_angular_momentum();
        void save(std::ofstream &f);
    public:
        void run();
        std::string get_name() {return integrator->name;}
        System(config &_c, bool _save_results = true);
        ~System() {delete integrator;}
};

template <typename Integrator>
System<Integrator>::System(config &_c, bool _save_results) : c{_c}, save_results{_save_results} {
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
    x = array(_data.data(), bodies * 3);
    v = array(_data.data() + bodies * 3, bodies * 3);
    m = array(_data.data() + bodies * 6, bodies);

    // banger ???
    if constexpr (std::is_base_of_v<FirstOrderODE, Integrator>) {
        integrator = new Integrator([&](fpoint_t t, array& Y, array& ret){this->compute_ac(t, Y, ret);}, data);
    } else if constexpr (std::is_base_of_v<SecondOrderODE, Integrator>) {
        integrator = new Integrator([&](fpoint_t t, array& Y, array& ret){this->compute_ac_order_2(t, Y, ret);}, [&](fpoint_t t, array& Y, array& ret){
            for (int i = 0; i < this->bodies * 3; i++) ret[i] = Y[i];
        }, x, v);
    } else {
       throw (std::bad_typeid());
    }
}

template <typename Integrator>
void System<Integrator>::compute_ac(fpoint_t t, array &Y, array &ret) {
    vecarray src_x = vecarray((vec3 *) Y.data_handle(), bodies);
    vecarray src_v = vecarray((vec3 *) Y.data_handle() + bodies, bodies);
    array src_m = array(Y.data_handle() + bodies * 6, bodies);
    vecarray ret_v = vecarray((vec3 *) ret.data_handle(), bodies);
    vecarray ret_a = vecarray((vec3 *) ret.data_handle() + bodies, bodies);
    array ret_m = array(ret.data_handle() + bodies * 6, bodies);

    for (int i = 0; i < bodies; i++) {
        vec3 a_i = vec3{0, 0, 0};

        for (int j = 0; j < i; j++) {
            vec3 r = src_x[j] - src_x[i];
            fpoint_t r_norm = r.norm();
            fpoint_t r3 = r_norm * r_norm * r_norm + eps_r; 
            a_i = a_i + r * gamma * src_m[i] * src_m[j] / r3;
        }

        for (int j = i+1; j < bodies; j++) {
            vec3 r = src_x[j] - src_x[i];
            fpoint_t r_norm = r.norm();
            fpoint_t r3 = r_norm * r_norm * r_norm + eps_r; 
            a_i = a_i + r * gamma * src_m[i] * src_m[j] / r3;
        }

        ret_a[i] = a_i / src_m[i];
        ret_v[i] = src_v[i];
        ret_m[i] = 0;
    } 
}

template <typename Integrator>
void System<Integrator>::compute_ac_order_2(fpoint_t t, array &_x, array &ret) {
    vecarray x_v = vecarray((vec3 *) _x.data_handle(), bodies);
    vecarray a_v = vecarray((vec3 *) ret.data_handle(), bodies);

    for (int i = 0; i < bodies; i++) {
        vec3 a_i = vec3{0, 0, 0};

        for (int j = 0; j < i; j++) {
            vec3 r = x_v[j] - x_v[i];
            fpoint_t r_norm = r.norm();
            fpoint_t r3 = r_norm * r_norm * r_norm + eps_r; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        for (int j = i+1; j < bodies; j++) {
            vec3 r = x_v[j] - x_v[i];
            fpoint_t r_norm = r.norm();
            fpoint_t r3 = r_norm * r_norm * r_norm + eps_r;
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        a_v[i] = a_i / m[i];
    }
}

template <typename Integrator>
void System<Integrator>::run() {
    std::ofstream f;
    f.open(integrator->name + ".txt", std::ios::out);
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
    f << "A_norm;K;U;K+U" << std::endl;;

    time_point<high_resolution_clock> t1 = high_resolution_clock::now();

    if (save_results) {
        compute_energies();
        compute_angular_momentum();
        save(f);
    } 
    while (integrator->t < c.target_t) {
        integrator->t += integrator->dt;
        integrator->step();
        integrator->steps++;
        if (save_results) {
            compute_energies();
            compute_angular_momentum();
            save(f);
        }
    }

    f.close();

    time_point<high_resolution_clock> t2 = high_resolution_clock::now();
    std::cout << integrator->name << ";\t";
    std::cout << (t2 - t1).count() / 1e9 << "s;\t"
              << integrator->steps << ";\t"
              << integrator->f_evals
              << std::endl;
}

template <typename Integrator>
void System<Integrator>::compute_energies() {
    size_t size = bodies * 3;
    vecarray x = vecarray((vec3 *) data.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) data.data_handle() + bodies, bodies);
    array m = array(data.data_handle() + bodies * 6, bodies);

    K = 0.0;
    U = 0.0;
    for (int i = 0; i < bodies; i++) {
        K += v[i].norm() * v[i].norm() * m[i] / 2.0;
        for (int j = i+1; j < m.size(); j++) {
            U += -gamma * m[i] * m[j] / (x[j] - x[i]).norm();
        }
    }
}

template <typename Integrator>
void System<Integrator>::compute_angular_momentum() {
    vecarray x = vecarray((vec3 *) data.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) data.data_handle() + bodies, bodies);
    A = vec3{0, 0, 0}; 
    for (int i = 0; i < bodies; i++) {
        A = A + v[i].prod(x[i]);
    }
}

template <typename Integrator>
void System<Integrator>::save(std::ofstream &f) {
    vecarray x = vecarray((vec3 *) data.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) data.data_handle() + bodies, bodies);
    array m = array(data.data_handle() + bodies * 6, bodies);

    f << std::setprecision(10)
      << bodies << ";" << integrator->steps << ";" << integrator->t << ";";
    for (int i = 0; i < bodies; i++) {
        f << m[i] << ";" 
          << x[i].x << ";" << x[i].y << ";" << x[i].z << ";"
          << v[i].x << ";" << v[i].y << ";" << v[i].z << ";";
    }

    f << A.norm() << ";" << K << ";" << U << ";" << K+U
      << std::endl;
}