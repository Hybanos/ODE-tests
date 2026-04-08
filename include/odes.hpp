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
    public:
        double t = 0;
        double dt = 0.001;
        int steps = 0;
        std::string name = "ODE";

        virtual void step() = 0;
        ODE(std::string _name) : name{_name} {}
};

class FirstOrderODE : public ODE {
    protected:
        ftype f;
        array Y;

        std::vector<double> _tmp;
        array tmp;
    public:
        FirstOrderODE(ftype f, array &Y0, std::string name) : ODE(name), f{f}, Y{Y0} {
            nd = Y.extent(0);
            _tmp.resize(nd);
            tmp = array(_tmp.data(), nd);
        }
};

class SecondOrderODE : public ODE {
    protected:
        ftype f;
        ftype g;

        array x;
        array v;
    public:
        SecondOrderODE(ftype _f, ftype _g, array &_x, array &_v, std::string name) : 
            ODE(name), f{_f}, g{_g}, x{_x}, v{_v} {
            nd = x.extent(0);
        }
};

class Euler : public FirstOrderODE {
    protected:
    public:
        void step();
        Euler(ftype f, array &Y0) : FirstOrderODE(f, Y0, "Euler") {}
};

class EulerSymplectic : public SecondOrderODE { 
    private:
        std::vector<double> _tmp1;
        std::vector<double> _tmp2;
        array tmp1;
        array tmp2;
    public:
        EulerSymplectic(ftype f, ftype g, array &x, array &v) : SecondOrderODE(f, g, x, v, "EulerSymplectic") {
            _tmp1.resize(nd);
            _tmp2.resize(nd);
            tmp1 = array(_tmp1.data(), nd);
            tmp2 = array(_tmp2.data(), nd);
        }
        void step();
};

class RK4 : public FirstOrderODE {
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
        RK4(ftype f, array &Y0) : FirstOrderODE(f, Y0, "RK4") {
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