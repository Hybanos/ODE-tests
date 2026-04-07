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
        double t = 0;
        double dt = 0.001;
        int steps = 0;
        int nd;
        ftype f;
        array Y;
    public:
        virtual void step() = 0;
        ODE(int _nd, ftype _f, array _Y0): nd{_nd}, f{_f}, Y{_Y0} {}
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
        RK4(int nd, ftype f, array X0) : ODE(nd, f, X0) {
            _k1.resize(nd);
            _k2.resize(nd);
            _k3.resize(nd);
            _k4.resize(nd);

            _tmp.resize(nd);

            k1 = array(_k1.data(), _k1.size());
            k2 = array(_k2.data(), _k2.size());
            k3 = array(_k3.data(), _k3.size());
            k4 = array(_k4.data(), _k4.size());

            tmp = array(_tmp.data(), _tmp.size());
        };
};