#include "odes.hpp"

void Euler::step() {
    f(t, Y, tmp);
    for (int i = 0; i < nd; i++) {
        Y[i] = Y[i] + tmp[i] * dt;
    }
    f_evals += 1;
}

void EulerSymplectic::step() {
    f(t, x, tmp1);
    for (int i = 0; i < nd; i++) {
        v[i] = v[i] + tmp1[i] * dt;
    }
    g(t, v, tmp2);
    for (int i = 0; i < nd; i++) {
        x[i] = x[i] + tmp2[i] * dt;
    } 
    // 2 ?
    f_evals += 1;
}

void LeapFrog::step() {
    f(t, x, tmp1);
    for (int i = 0; i < nd; i++) {
        x[i] = x[i] + v[i] * dt + tmp1[i] * 0.5 * dt * dt;
    }

    f(t, x, tmp2);
    for (int i = 0; i < nd; i++) {
        v[i] = v[i] + (tmp1[i] + tmp2[i]) * dt * 0.5;
    }
    f_evals += 2;
}

void RK2::step() {
    f(t, Y, k1);
    for (int i = 0; i < nd; i++) {
        tmp[i] = Y[i] + k1[i] * theta * dt;
    }

    f(t + theta * dt, tmp, k2);
    for (int i = 0; i < nd; i++) {
        Y[i] = Y[i] + ((1 - 0.5 * theta) * k1[i] + 0.5 * theta * k2[i]) * dt;
    }
    f_evals += 2;
}

void RK4::step() {
    f(t, Y, k1);
    for (int i = 0; i < nd; i++) {
        tmp[i] = Y[i] + k1[i] * 0.5 * dt;
    }

    f(t + dt/2, tmp, k2);
    for (int i = 0; i < nd; i++) {
        tmp[i] = Y[i] + k2[i] * 0.5 * dt;
    }

    f(t + dt/2, tmp, k3);
    for (int i = 0; i < nd; i++) {
        tmp[i] = Y[i] + k3[i] * dt;
    }

    f(t + dt, tmp, k4);

    for (int i = 0; i < nd; i++) {
        Y[i] = Y[i] + (k1[i] + k2[i]*2.0 + k3[i]*2.0 + k4[i]) * dt / 6.0;
    }
    f_evals += 4;
}

void RK45::step() {
    fpoint_t sum_e;
    do {
        f(t, Y, k1);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + k1[i] * dt * B[1][0];
        }

        f(t, tmp, k2);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * B[2][0] + k2[i] * B[2][1]) * dt;
        }

        f(t, tmp, k3);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * B[3][0] + k2[i] * B[3][1] + k3[i] * B[3][2]) * dt;
        }

        f(t, tmp, k4);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * B[4][0] + k2[i] * B[4][1] + k3[i] * B[4][2] + k4[i] * B[4][3]) * dt;
        }

        f(t, tmp, k5);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * B[5][0] + k2[i] * B[5][1] + k3[i] * B[5][2] + k4[i] * B[5][3] + k5[i] * B[5][4]) * dt;
        }

        f(t, tmp, k6);

        sum_e = 0;
        for (int i = 0; i < nd; i++) {
            sum_e += dt * std::abs(k1[i] * (ch[0] - c[0]) + k2[i] * (ch[1] - c[1]) + k3[i] * (ch[2] - c[2]) + k4[i] * (ch[3] - c[3]) + k5[i] * (ch[4] - c[4]) + k6[i] * (ch[5] - c[5]));
        }
        sum_e = std::abs(sum_e) * nd;

        dt = 0.9 * dt * std::pow(eps / sum_e, 1.0/5.0);
        f_evals += 6;
    } while (sum_e > eps);

    for (int i = 0; i < nd; i++) {
        Y[i] = Y[i] + (k1[i] * ch[0] + k2[i] * ch[1] + k3[i] * ch[2] + k4[i] * ch[3] + k5[i] * ch[4] + k6[i] * ch[5]) * dt;
    }

}

void DOP853::step() {
    using namespace dop853_coefs;

    // NOTE: ref impl takes an educated guess on the value of dt for step 0,
    // we're skipping that i'm lazy

    fpoint_t err = 0;
    do {
        f(t, Y, k1);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a21) * dt;
        }

        f(t, tmp, k2);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a31 + k2[i] * a32) * dt;
        }

        f(t, tmp, k3);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a41 + k3[i] * a43) * dt;
        }

        f(t, tmp, k4);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a51 + k3[i] * a53 + k4[i] * a54) * dt;
        }

        f(t, tmp, k5);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a61 + k4[i] * a64 + k5[i] * a65) * dt;
        }

        f(t, tmp, k6);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a71 + k4[i] * a74 + k5[i] * a75 + k6[i] * a76) * dt;
        }

        f(t, tmp, k7);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a81 + k4[i] * a84 + k5[i] * a85 + k6[i] * a86 + k7[i] * a87) * dt;
        }

        f(t, tmp, k8);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a91 + k4[i] * a94 + k5[i] * a95 + k6[i] * a96 + k7[i] * a97
                    + k8[i] * a98) * dt;
        }

        f(t, tmp, k9);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a101 + k4[i] * a104 + k5[i] * a105 + k6[i] * a106 
                    + k7[i] * a107 + k8[i] * a108 + k9[i] * a109) * dt;
        }

        f(t, tmp, k10);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a111 + k4[i] * a114 + k5[i] * a115 + k6[i] * a116 
                    + k7[i] * a117 + k8[i] * a118 + k9[i] * a119 + k10[i] * a1110) * dt;
        }

        // k2 reuse !! not a bug !! 
        f(t, tmp, k2);
        for (int i = 0; i < nd; i++) {
            tmp[i] = Y[i] + (k1[i] * a121 + k4[i] * a124 + k5[i] * a125 + k6[i] * a126 
                    + k7[i] * a127 + k8[i] * a128 + k9[i] * a129 + k10[i] * a1210 + k2[i] * a1211) * dt;
        }

        // same for k3 !!1!
        f(t, tmp, k3);
        for (int i = 0; i < nd; i++) {
            k4[i] = k1[i] * b1 + k6[i] * b6 + k7[i] * b7 + k8[i] * b8 + k9[i] * b9 + k10[i] * b10 + k2[i] * b11 + k3[i] * b12;
            k5[i] = Y[i] + k4[i] * dt;
        } 

        // error estimation
        err = 0.0;
        fpoint_t err_2 = 0.0;
        for (int i = 0; i < nd; i++) {
            fpoint_t sk = a_tol + r_tol * std::max(std::abs(Y[i]), std::abs(k5[i]));
            fpoint_t err_i = k4[i] - k1[i] * bhh1 - k9[i] * bhh2 - k3[i] * bhh3;
            err_2 = err_2 + (err_i / sk) * (err_i / sk);
            err_i = k1[i] * er1 + k6[i] * er6 + k7[i] * er7 + k8[i] * er8 + k9[i] * er9 + k10[i] * er10 + k2[i] * er11 + k3[i] * er12;
            err = err + (err_i / sk) * (err_i / sk);
        }

        fpoint_t deno = err + 0.01 * err_2;
        if (deno <= 0) deno = 1.0;
        err = dt * err * std::sqrt(1.0 / (nd * deno));
        // new dt 
        fpoint_t fac11 = std::pow(err, 1.0 / 8.0 - beta * 0.2);
        fpoint_t fac = fac11 / std::pow(facold, beta);
        fac = std::max((fpoint_t) 1.0/6, std::min(facc1, fac/safe));
        // std::cout << steps << " " << dt << " " << err << std::endl;

        f_evals += 12;
        if (err < 1.0) {
            facold = std::max(err, (fpoint_t) 1e-4);
            dt = dt / fac;
            break;
        }
        dt = dt / std::min(facc1, fac11 / safe);
    } while (1);

    // we skip both stiffness detection and dense output prep
    for (int i = 0; i < nd; i++) {
        Y[i] = k5[i];
    }
}

ref_ftype haha;
std::ofstream fout;
int ndims;

void dop_ref_solout(int *nr, double *xold, double *X, double *Y, 
                    int *n, double *con, int *icomp, int *nrd, 
                    int *rpar, int *ipar, int *irtrn, double *xout) {

    int bodies = ndims / 7;
    vecarray x = vecarray((vec3 *) Y, bodies);
    vecarray v = vecarray((vec3 *) Y + bodies, bodies);
    array m = array(Y + bodies * 6, bodies);

    fout << std::setprecision(10)
      << bodies << ";" << *nr << ";" << *X << ";";
    for (int i = 0; i < bodies; i++) {
        fout << m[i] << ";" 
          << x[i].x << ";" << x[i].y << ";" << x[i].z << ";"
          << v[i].x << ";" << v[i].y << ";" << v[i].z << ";";
    }

    fout << -1 << ";" << -1 << ";" << -1 << ";" << -1 
      << std::endl;
}

void dop_ref_f_wrapper(int *n, double *t, double * X, double * F, double *idk, int *idk2) {
    haha(n, t, X, F, idk, idk2);
}

void DOP853_ref::step() {
    int n = Y.extent(0);
    ndims = n;
    haha = [&](int *n, double *t, double * X, double * F, double *idk, int *idk2){
        array x(X, *n);
        array f(F, *n);
        this->f(this->t, x, f);
    };
    double x = 0.0;
    double *y = Y.data_handle();
    double xend = 10.0;
    double rtol = R_TOL;
    double atol = A_TOL;
    int itol = 0;
    int iout = 1;
    int lwork = 11 * n + 21;
    double *work = new double[lwork];
    int liwork = 21;
    int *iwork = new int[liwork];
    double rpar = 0;
    int ipar = 0;
    int idid = 0;

    for (int i = 0; i< lwork; i++) work[i] = 0.0;
    for (int i = 0; i< liwork; i++) iwork[i] = 0;

    // std::cout << (void *) dop_ref_f_wrapper << std::endl;

    fout.open("DOP853-ref.txt", std::ios::app);
    dop853(&n, (void *)dop_ref_f_wrapper, &x, Y.data_handle(), 
        &xend, &rtol, &atol, &itol, (void *)dop_ref_solout, 
        &iout, work, &lwork, iwork, &liwork, 
        &rpar, &ipar, &idid);
    t = 10;
    fout.close();
}
