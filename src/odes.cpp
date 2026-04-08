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