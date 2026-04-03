#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <vector>

#include <mdspan/mdspan.hpp>

#include "vec3.hpp"

namespace stdex = Kokkos;

using extents = stdex::extents<int, stdex::dynamic_extent>;
using array = stdex::mdspan<double, extents>;
using vecarray = stdex::mdspan<vec3, extents>;

class System {
    protected:
        double gamma = 1;

        double t = 0.0;
        double dt = 0.05;
        double target_t;
        int steps = 0;

        std::vector<double> m;
        std::vector<vec3> x;
        std::vector<vec3> v;

        // acceleration buffer
        std::vector<vec3> a;

        std::vector<vec3> _data;
        array data;

        // kinetic and potential energies
        double K, U;

        void compute_acc_md(double t, array &X, array &ret);
        void compute_accelerations(std::vector<vec3> &a, std::vector<vec3> &x);
        void compute_energies();
        virtual void step() = 0;
        void save(std::ofstream &f);
    public:
        std::string name;

        System(std::string _name, double _target_t, int bodies, int seed);
        void run();
};