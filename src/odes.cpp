#include "odes.hpp"

void Euler::step() {
    f(t, Y, tmp);
    for (int i = 0; i < nd; i++) {
        Y[i] = Y[i] + tmp[i] * dt;
    }
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
        Y[i] = Y[i] + (k1[i] + k2[i]*2 + k3[i]*2 + k4[i]) * dt / 6;
    }
}

void RK45::step() {
    double sum_e = 0;
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
            sum_e += k1[i] * (ch[0] - c[0]) + k2[i] * (ch[1] - c[1]) + k3[i] * (ch[2] - c[2]) + k4[i] * (ch[3] - c[3]) + k5[i] * (ch[4] - c[4]) + k6[i] * (ch[5] - c[5]);
        }
        sum_e = std::abs(sum_e) / nd;

        dt = 0.9 * dt * std::pow(eps / sum_e, 1.0/5.0);
    std::cout << dt << std::endl;
    } while (sum_e > eps);

    for (int i = 0; i < nd; i++) {
        Y[i] = Y[i] + (k1[i] * ch[0] + k2[i] * ch[1] + k3[i] * ch[2] + k4[i] * ch[3] + k5[i] * ch[4] + k6[i] * ch[5]) * dt;
    }

    if (dt > base_dt) dt = base_dt * 1;
}