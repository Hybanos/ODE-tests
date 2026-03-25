#pragma once

#include <cmath>

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
        Exact(int steps);
};

class Euler : public System {
    public:
        void step();
        Euler(int steps) : System{"Euler", steps} {}
};

class RK2 : public System {
    private:
        double theta;
    public:
        void step();
        RK2(int steps, double _theta) : System{std::string("RK2") + std::to_string((int) (1 / _theta)), steps}, theta{_theta} {}
};