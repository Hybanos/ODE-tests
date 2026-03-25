#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

#include "vec3.hpp"

class System {
    protected:
        int max_t;

        double gamma = 1;

        double t = 0.0;
        double dt = 0.001;
        int steps = 0;

        double m1 = 2.0;
        double m2 = 1.0;

        vec3 p1 = {-1, 0, 0};       
        vec3 p2 = { 1, 0, 0};       

        vec3 v1 = {0,-0.6 / m1, 0};
        vec3 v2 = {0, 0.6, 0};

        virtual void step() = 0;
        void print();
        void save(std::ofstream &f);
    public:
        std::string name;

        System(std::string _name, int _max_t) : name{_name}, max_t{_max_t} {};
        void run();
};