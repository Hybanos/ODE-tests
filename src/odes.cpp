#include "odes.hpp"

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