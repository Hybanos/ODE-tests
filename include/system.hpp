#include <iostream>
#include <cmath>

#include "vec3.hpp"

class System {
    private:
        int max_t;

        double gamma = 1;

        double t = 0.0;
        double dt = 1.0;

        double m1 = 1.0;
        double m2 = 1.0;

        vec3 p1 = {-100, 0, 0};       
        vec3 p2 = { 100, 0, 0};       

        vec3 v1 = {0, -1, 0};
        vec3 v2 = {0,  1, 0};

        void step();
        void print();
    public:
        System(int _max_t) : max_t{_max_t} {};
        void run();
};