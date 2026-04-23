#pragma once

#include <cmath>
#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>

#include <mdspan/mdspan.hpp>

#include "vec3.hpp"
#include "dop853coefs.hpp"
#include "defs.hpp"

namespace stdex = Kokkos;

using extents = stdex::extents<int, stdex::dynamic_extent>;
using array = stdex::mdspan<fpoint_t, extents>;
using vecarray = stdex::mdspan<vec3, extents>;
using ftype = std::function<void(fpoint_t t, array& Y, array& ret)>;

using ref_ftype = std::function<void(int *n, double *t, double * X, double * F, double *idk, int *idk2)>;
using ref_solout_ftype = std::function<void(int *nr, double *xold, double *x, double *Y, 
                    int *n, double *con, int *icomp, int *nd, 
                    int *rpar, int *ipar, int *irtrn, double *xout)>;

extern "C" {
    void dop853(int *n, void * f,
        double *x, double *Y, double *xend, double *rtol, double *atol, int *itol,
        void * solout,
        int *iout, double *work, int *lwork, int *iwork, int *liwork, double *rpar, int *ipar, int *idid
    );
}

class ODE {
    protected:
        int nd;
    public:
        fpoint_t t = 0;
        fpoint_t dt = 0.001;
        int steps = 0;
        int f_evals = 0;
        std::string name = "ODE";


        virtual void step() = 0;
        ODE(std::string _name) : name{_name} {}
};

class FirstOrderODE : public ODE {
    protected:
        ftype f;
        array Y;

        std::vector<fpoint_t> _tmp;
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
        std::vector<fpoint_t> _tmp1;
        std::vector<fpoint_t> _tmp2;
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

class LeapFrog : public SecondOrderODE { 
    private:
        std::vector<fpoint_t> _tmp1;
        std::vector<fpoint_t> _tmp2;
        array tmp1;
        array tmp2;
    public:
        LeapFrog(ftype f, ftype g, array &x, array &v) : SecondOrderODE(f, g, x, v, "LeapFrog") {
            _tmp1.resize(nd);
            _tmp2.resize(nd);
            tmp1 = array(_tmp1.data(), nd);
            tmp2 = array(_tmp2.data(), nd);
        }
        void step();
};

class RK2 : public FirstOrderODE {
    private:
        fpoint_t theta = 0.5;
        
        std::vector<fpoint_t> _k1;
        std::vector<fpoint_t> _k2;
        std::vector<fpoint_t> _tmp;

        array k1;
        array k2;
        array tmp;
    public:
        void step();
        RK2(ftype f, array &Y0) : FirstOrderODE(f, Y0, std::string("RK2") + std::to_string(2)) {
            _k1.resize(nd);
            _k2.resize(nd);
            _tmp.resize(nd);

            k1 = array(_k1.data(), nd);
            k2 = array(_k2.data(), nd);
            tmp = array(_tmp.data(), nd);
        }
};

class RK4 : public FirstOrderODE {
    private:
        std::vector<fpoint_t> _k1;
        std::vector<fpoint_t> _k2;
        std::vector<fpoint_t> _k3;
        std::vector<fpoint_t> _k4;

        std::vector<fpoint_t> _tmp;

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

class RK45 : public FirstOrderODE {
    private:
        fpoint_t base_dt;
        fpoint_t eps = A_TOL;
        fpoint_t A[6] = {2.0/9.0, 1.0/3.0, 3.0/4.0, 1, 5.0/6.0};
        fpoint_t B[6][5] = {
            {0,           0,           0,           0,         0},
            {2.0/9.0,     0,           0,           0,         0},
            {1.0/12.0,    1.0/4.0,     0,           0,         0},
            {69.0/128.0, -243.0/128.0, 135.0/64.0,  0,         0},
            {-17.0/12.0,  27.0/4.0,   -27.0/5.0,    16.0/15.0, 0},
            {65.0/432.0, -5.0/16.0,    13.0/16.0,   4.0/27.0,  5.0/144.0}
        };
        fpoint_t c[6] = {1.0/9.0, 0, 9.0/20.0, 16.0/45.0, 1.0/12.0};
        fpoint_t ch[6] = {47.0/450.0, 0, 12.0/25.0, 32.0/225.0, 1.0/30.0, 6.0/25.0};

        std::vector<fpoint_t> _k1;
        std::vector<fpoint_t> _k2;
        std::vector<fpoint_t> _k3;
        std::vector<fpoint_t> _k4;
        std::vector<fpoint_t> _k5;
        std::vector<fpoint_t> _k6;
        std::vector<fpoint_t> _tmp;

        array k1;
        array k2;
        array k3;
        array k4;
        array k5;
        array k6;
        array tmp;

    public:
        void step();
        RK45(ftype f, array &Y0) : FirstOrderODE(f, Y0, "RK45") {
            base_dt = dt;

            _k1.resize(nd);
            _k2.resize(nd);
            _k3.resize(nd);
            _k4.resize(nd);
            _k5.resize(nd);
            _k6.resize(nd);
            _tmp.resize(nd);

            k1 = array(_k1.data(), nd);
            k2 = array(_k2.data(), nd);
            k3 = array(_k3.data(), nd);
            k4 = array(_k4.data(), nd);
            k5 = array(_k5.data(), nd);
            k6 = array(_k6.data(), nd);
            tmp = array(_tmp.data(), nd);
        }
};

class DOP853 : public FirstOrderODE {
    private:
        fpoint_t base_dt;    
        fpoint_t a_tol = A_TOL;
        fpoint_t r_tol = R_TOL;

        fpoint_t beta = 0;

        fpoint_t safe = 0.9;
        fpoint_t fac1 = 0.3333;
        fpoint_t facc1 = 1.0 / fac1;
        fpoint_t facold = 1e-4;

        std::vector<fpoint_t>  _k1;
        std::vector<fpoint_t>  _k2;
        std::vector<fpoint_t>  _k3;
        std::vector<fpoint_t>  _k4;
        std::vector<fpoint_t>  _k5;
        std::vector<fpoint_t>  _k6;
        std::vector<fpoint_t>  _k7;
        std::vector<fpoint_t>  _k8;
        std::vector<fpoint_t>  _k9;
        std::vector<fpoint_t> _k10;
        std::vector<fpoint_t> _tmp;

        array  k1;
        array  k2;
        array  k3;
        array  k4;
        array  k5;
        array  k6;
        array  k7;
        array  k8;
        array  k9;
        array k10;
        array tmp;
    public:
        void step();
        DOP853(ftype f, array &Y0) : FirstOrderODE(f, Y0, "DOP853") {
            base_dt = dt;

            _k1.resize(nd);
            _k2.resize(nd);
            _k3.resize(nd);
            _k4.resize(nd);
            _k5.resize(nd);
            _k6.resize(nd);
            _k7.resize(nd);
            _k8.resize(nd);
            _k9.resize(nd);
            _k10.resize(nd);
            _tmp.resize(nd);

            k1 = array(_k1.data(), nd);
            k2 = array(_k2.data(), nd);
            k3 = array(_k3.data(), nd);
            k4 = array(_k4.data(), nd);
            k5 = array(_k5.data(), nd);
            k6 = array(_k6.data(), nd);
            k7 = array(_k7.data(), nd);
            k8 = array(_k8.data(), nd);
            k9 = array(_k9.data(), nd);
            k10 =array(_k10.data(), nd);
            tmp =array(_tmp.data(), nd);
        }
        
};

class DOP853_ref : public FirstOrderODE {
    private:
    public:
        bool save_results = false;
        void step();
        DOP853_ref(ftype f, array &Y0) : FirstOrderODE(f, Y0, "DOP853-ref") {}
};