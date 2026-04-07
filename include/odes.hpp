#pragma once

#include <cmath>
#include <vector>
#include <functional>

#include <mdspan/mdspan.hpp>

#include "vec3.hpp"
#include "system.hpp"
#include "dop853coefs.hpp"

using ftype = std::function<void(double t, array& Y, array& ret)>;

class ODE {
    protected:
        int nd;
        ftype f;
        array Y;

        std::vector<double> _tmp;
        array tmp;
    public:
        double t = 0;
        double dt = 0.001;
        int steps = 0;
        std::string name = "ODE";

        virtual void step() = 0;
        ODE(ftype _f, array &_Y0, std::string _name) : f{_f}, Y{_Y0}, name{_name} {
            nd = Y.extent(0);
            _tmp.resize(nd);
            tmp = array(_tmp.data(), nd);
        }
};

class Euler : public ODE {
    public:
        void step();
        Euler(ftype f, array &Y0) : ODE(f, Y0, "Euler") {}
};

class RK4 : public ODE {
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
        RK4(ftype f, array &Y0) : ODE(f, Y0, "RK4") {
            _k1.resize(nd);
            _k2.resize(nd);
            _k3.resize(nd);
            _k4.resize(nd);

            _tmp.resize(nd);

            k1 = array(_k1.data(), nd);
            k2 = array(_k2.data(), nd);
            k3 = array(_k3.data(), nd);
            k4 = array(_k4.data(), nd);

            tmp = array(_tmp.data(), nd);
        };
};