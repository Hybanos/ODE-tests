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

template <typename Integrator>
class System {
    private:
        // initial configuration
        config c;

        Integrator *integrator = nullptr;

        int bodies;
        double gamma = 1;
        // double target_t = 0.001;
        double target_t = 90;

        std::vector<double> _data;
        array data;
        array x;
        array v;
        array m;

        // kinetic and potential energies
        double K, U;

        void compute_ac(double t, array &Y, array &ret);
        void compute_ac_order_2(double t, array &x, array &ret);
        void compute_energies();
        void save(std::ofstream &f);
    public:
        void run();
        std::string get_name() {return integrator->name;}
        System(config &_c);
        ~System() {delete integrator;}
};

template <typename Integrator>
System<Integrator>::System(config &_c) : c{_c} {
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
        integrator = new Integrator([&](double t, array& Y, array& ret){this->compute_ac(t, Y, ret);}, data);
    } else if constexpr (std::is_base_of_v<SecondOrderODE, Integrator>) {
        integrator = new Integrator([&](double t, array& Y, array& ret){this->compute_ac_order_2(t, Y, ret);}, [&](double t, array& Y, array& ret){
            for (int i = 0; i < this->bodies * 3; i++) ret[i] = Y[i];
        }, x, v);
    } else {
       throw (std::bad_typeid());
    }
}

template <typename Integrator>
void System<Integrator>::compute_ac(double t, array &Y, array &ret) {
    vecarray x_v = vecarray((vec3 *) Y.data_handle(), bodies);
    vecarray v_v = vecarray((vec3 *) Y.data_handle() + bodies, bodies);
    vecarray ret_v = vecarray((vec3 *) ret.data_handle(), bodies);
    vecarray ret_a = vecarray((vec3 *) ret.data_handle() + bodies, bodies);

    for (int i = 0; i < bodies; i++) {
        vec3 a_i = vec3{0, 0, 0};

        for (int j = 0; j < i; j++) {
            vec3 r = x_v[j] - x_v[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        for (int j = i+1; j < bodies; j++) {
            vec3 r = x_v[j] - x_v[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        ret_a[i] = a_i / m[i];
        ret_v[i] = v_v[i];
    } 
}

template <typename Integrator>
void System<Integrator>::compute_ac_order_2(double t, array &_x, array &ret) {
    vecarray x_v = vecarray((vec3 *) _x.data_handle(), bodies);
    vecarray a_v = vecarray((vec3 *) ret.data_handle(), bodies);

    for (int i = 0; i < bodies; i++) {
        vec3 a_i = vec3{0, 0, 0};

        for (int j = 0; j < i; j++) {
            vec3 r = x_v[j] - x_v[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
            a_i = a_i + r * gamma * m[i] * m[j] / r3;
        }

        for (int j = i+1; j < bodies; j++) {
            vec3 r = x_v[j] - x_v[i];
            double r_norm = r.norm();
            double r3 = r_norm * r_norm * r_norm; 
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
    f << "K;U;K+U" << std::endl;;

    time_point<high_resolution_clock> t1 = high_resolution_clock::now();

    compute_energies();
    save(f);
    while (integrator->t < target_t) {
        integrator->t += integrator->dt;
        integrator->step();
        compute_energies();
        save(f);
    }

    f.close();

    time_point<high_resolution_clock> t2 = high_resolution_clock::now();
    std::cout << integrator->name << ": ";
    std::cout << (t2 - t1).count() / 1e9 << "s" << std::endl;
}

template <typename Integrator>
void System<Integrator>::compute_energies() {
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

template <typename Integrator>
void System<Integrator>::save(std::ofstream &f) {
    vecarray x = vecarray((vec3 *) data.data_handle(), bodies);
    vecarray v = vecarray((vec3 *) data.data_handle() + bodies, bodies);
    array m = array(data.data_handle() + bodies * 6, bodies);

    f << bodies << ";" << integrator->steps << ";" << integrator->t << ";";
    for (int i = 0; i < bodies; i++) {
        f << m[i] << ";" 
          << x[i].x << ";" << x[i].y << ";" << x[i].z << ";"
          << v[i].x << ";" << v[i].y << ";" << v[i].z << ";";
    }

    f << K << ";" << U << ";" << K+U
      << std::endl;
}