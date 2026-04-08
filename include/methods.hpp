#pragma once

#include <cmath>
#include <vector>

#include <mdspan/mdspan.hpp>

#include "vec3.hpp"
#include "system.hpp"
#include "dop853coefs.hpp"

class Exact : public System {
    private:
        vec3 barycenter_pos; 
        vec3 barycenter_speed; 

        vec3 h;

        vec3 vec_eccentricity;

        double nu;
        double semi_major_axis;
        double eccentricity;
        double E_0;
        double mean_motion;

        double solve_true_anomaly();
    public:
        void step();
        Exact(double target_t, int bodies, int seed);
};

class Euler : public System {
    public:
        void step();
        Euler(double target_t, int bodies, int seed=0) : System{"Euler", target_t, bodies, seed} {}
};

class EulerSwapped : public System {
    public:
        void step();
        EulerSwapped(double target_t, int bodies, int seed=0) : System{"EulerSwapped", target_t, bodies, seed} {}
};

class Leapfrog : public System {
    public:
        void step();
        Leapfrog(double target_t, int bodies, int seed=0) : System{"Leapfrog", target_t, bodies, seed} {}
};

class RK2 : public System {
    private:
        double theta;
    public:
        void step();
        RK2(double target_t, int bodies, int seed=0, double _theta=1) : System{std::string("RK2") + std::to_string((int) (1 / _theta)), target_t, bodies, seed}, theta{_theta} {}
};

class RK4 : public System {
    public:
        void step();
        RK4(double target_t, int bodies, int seed=0) : System{"RK4", target_t, bodies, seed} {}
};

class RK4_md : public System {
    private:
        std::vector<double> _k1;
        std::vector<double> _k2;
        std::vector<double> _k3;
        std::vector<double> _k4;

        std::vector<double> _tmp;

        array k1;
        array k2;
        array k3;
        array k4;

        array tmp;
    public:
        void step();
        RK4_md(double target_t, int bodies, int seed=0) : System{"RK4-md", target_t, bodies, seed} {
            _k1.resize(data.extent(0));
            _k2.resize(data.extent(0));
            _k3.resize(data.extent(0));
            _k4.resize(data.extent(0));

            _tmp.resize(data.extent(0));

            k1 = array(_k1.data(), _k1.size());
            k2 = array(_k2.data(), _k2.size());
            k3 = array(_k3.data(), _k3.size());
            k4 = array(_k4.data(), _k4.size());

            tmp = array(_tmp.data(), _tmp.size());
        }
};

class RK45 : public System {
    private:
        double base_dt;
        double eps = 1e-3;
        double A[6] = {2.0/9.0, 1.0/3.0, 3.0/4.0, 1, 5.0/6.0};
        double B[6][5] = {
            {0,           0,           0,           0,         0},
            {2.0/9.0,     0,           0,           0,         0},
            {1.0/12.0,    1.0/4.0,     0,           0,         0},
            {69.0/128.0, -243.0/128.0, 135.0/64.0,  0,         0},
            {-17.0/12.0,  27.0/4.0,   -27.0/5.0,    16.0/15.0, 0},
            {65.0/432.0, -5.0/16.0,    13.0/16.0,   4.0/27.0,  5.0/144.0}
        };
        double c[6] = {1.0/9.0, 0, 9.0/20.0, 16.0/45.0, 1.0/12.0};
        double ch[6] = {47.0/450.0, 0, 12.0/25.0, 32.0/225.0, 1.0/30.0, 6.0/25.0};
    public:
        void step();
        RK45(double target_t, int bodies, int seed=0) : System{"RK45", target_t, bodies, seed} {base_dt = dt;}
};

class DOP853 : public System {
    private:
        double base_dt;
        double a_tol = 1e-3;
        double r_tol = 1e-6;

        double beta = 0;

        double facold = 1e-4;
    public:
        void step();
        DOP853(double target_t, int bodies, int seed=0) : System{"DOP853", target_t, bodies, seed} {base_dt = dt;}
};

class LinearMultistep : public System {
    private:
        int back_steps;
        std::vector<vec3> prev;
    public:
        void step();
        LinearMultistep(double target_t, int _back_steps, int bodies, int seed);
};