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

        // acceleration buffer
        std::vector<vec3> a;

        // kinetic and potential energies
        double K, U;

        vec3 compute_acceleration(vec3 x1, vec3 x2);
        void compute_accelerations(std::vector<vec3> &a, std::vector<vec3> &x);
        void compute_energies();
        virtual void step() = 0;
        void save(std::ofstream &f);
    public:
        std::string name;

        System(std::string _name, int _max_t, int bodies, int seed);
        void run();
};