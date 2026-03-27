#pragma once

#include <cmath>
#include <vector>

#include "vec3.hpp"
#include "system.hpp"

class Exact : public System {
    private:
        vec3 barycenter_pos; 
        vec3 barycenter_speed; 

        double semi_major_axis;
        double eccentricity;
        double E_0;
        double mean_movement;

        vec3 compute_barycenter(double t);
        vec3 compute_pos(double t);
    public:
        void step();
        Exact(int max_steps, int bodies, int seed);
};

class Euler : public System {
    public:
        void step();
        Euler(int max_steps, int bodies, int seed=0) : System{"Euler", max_steps, bodies, seed} {}
};

class EulerSwapped : public System {
    public:
        void step();
        EulerSwapped(int max_steps, int bodies, int seed=0) : System{"EulerSwapped", max_steps, bodies, seed} {}
};

class Leapfrog : public System {
    public:
        void step();
        Leapfrog(int max_steps, int bodies, int seed=0) : System{"Leapfrog", max_steps, bodies, seed} {}
};

class RK2 : public System {
    private:
        double theta;
    public:
        void step();
        RK2(int max_steps, int bodies, int seed=0, double _theta=1) : System{std::string("RK2") + std::to_string((int) (1 / _theta)), max_steps, bodies, seed}, theta{_theta} {}
};

class RK4 : public System {
    public:
        void step();
        RK4(int max_steps, int bodies, int seed=0) : System{"RK4", max_steps, bodies, seed} {}
};

class LinearMultistep : public System {
    private:
        int back_steps;
        vec3 prev1;
        vec3 prev2;
    public:
        void step();
        LinearMultistep(int max_steps, int _back_steps, int bodies, int seed);
};