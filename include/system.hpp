#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>

#include "vec3.hpp"

class System {
    protected:
        int max_t;

        double gamma = 1;

        double t = 0.0;
        double dt = 0.001;
        int steps = 0;

        std::vector<double> m;
        std::vector<vec3> x;
        std::vector<vec3> v;

        // double m1 = 2.0;
        // double m2 = 1.0;

        // vec3 p1 = {-1, 0, 0};       
        // vec3 p2 = { 1, 0, 0};       

        // vec3 v1 = {0,-0.6 / m1, 0};
        // vec3 v2 = {0, 0.6, 0};

        // kinetic and potential energies
        double K, U;

        vec3 compute_acceleration(vec3 x1, vec3 x2);
        void compute_energies();
        virtual void step() = 0;
        void save(std::ofstream &f);
    public:
        std::string name;

        System(std::string _name, int _max_t, int bodies, int seed);
        void run();
};