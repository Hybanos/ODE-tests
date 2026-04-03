#include "methods.hpp"

Exact::Exact(double target_t, int bodies, int seed=0) : System("Exact", target_t, bodies, seed) {

    if (m.size() > 2) {exit(1);}

    double &m1 = m[0];
    double &m2 = m[1];
    vec3 &x1 = x[0];
    vec3 &x2 = x[1];
    vec3 &v1 = v[0];
    vec3 &v2 = v[1];

    barycenter_pos = (x1 * m1 + x2 * m2) / (m1 + m2);
    barycenter_speed = (v1 * m1 + v2 * m2) / (m1 + m2);

    vec3 r = x1 - x2;
    vec3 v = v1 - v2;
    h = r.prod(v);

    double eps = (v.norm() * v.norm()) / 2 - gamma * (m1 + m2) / r.norm();
    std::cout << "epsilon: " << eps << std::endl;

    semi_major_axis = - gamma * (m1 + m2) / (2 * eps);
    std::cout << "semi major axis: " << semi_major_axis << std::endl;

    vec_eccentricity = v.prod(h) / (gamma * (m1 + m2)) - r / r.norm();
    eccentricity = vec_eccentricity.norm();
    std::cout << "eccentricity: " << eccentricity << std::endl;

    double cos_nu = vec_eccentricity.dot(r) / (eccentricity * r.norm());
    nu = std::acos(cos_nu);
    if (r.dot(v) < 0) nu = 2*M_PI - nu;
    std::cout << "initial true anomaly: " << nu << std::endl;

    mean_motion = std::sqrt(gamma * (m1 + m2) / (semi_major_axis * semi_major_axis * semi_major_axis));
    std::cout << "mean motion: " << mean_motion << std::endl;
}

double Exact::solve_true_anomaly() {
    double e = 1e-9;
    int max = 10000;

    // double tan_E_2 = std::sqrt((1 - eccentricity) / (1 + eccentricity)) * std::tan(nu / 2);
    // double E_0 = std::atan(tan_E_2 / 2);
    
    // // double E_0 = 2 * std::atan2(std::sqrt(1 - eccentricity) * std::sin(nu / 2), std::sqrt(1 + eccentricity) * std::cos(nu / 2));
    // double M_0 = E_0 - eccentricity * std::sin(E_0);

    // std::cout << E_0 << " " << M_0 << std::endl;

    // double M_0 = nu - 2 * eccentricity * sin(nu) +
    //              (0.75 * eccentricity * eccentricity + 1 / 8 * eccentricity * eccentricity * eccentricity * eccentricity) * sin(2 * nu) -
    //              1 / 3 * eccentricity * eccentricity * eccentricity * sin(3 * nu) +
    //              5 / 32 * eccentricity * eccentricity * eccentricity * eccentricity * sin(4 * nu);

    double tmp_mean_anomaly = mean_motion * t;
    double anm_ecc = eccentricity < 0.8 ? tmp_mean_anomaly : M_PI;
    double f = anm_ecc - eccentricity * std::sin(anm_ecc) - tmp_mean_anomaly;
    int ite = 0;

    while (std::abs(f) > e) {
        f = anm_ecc - eccentricity * sin(anm_ecc) - tmp_mean_anomaly;
        anm_ecc = anm_ecc - f / (1.0 - eccentricity * std::cos(anm_ecc));

        ite += 1;
        if (ite > max) {
            std::cout << "failed to converge" << std::endl;
            break;
        }
    }

    return std::fmod(anm_ecc + M_PI * 2, M_PI * 2);
}

void Exact::step() {
    double &m1 = m[0];
    double &m2 = m[1];
    vec3 &x1 = x[0];
    vec3 &x2 = x[1];
    vec3 &v1 = v[0];
    vec3 &v2 = v[1];

    double true_anomaly = solve_true_anomaly();

    double x = semi_major_axis * (cos(true_anomaly) - eccentricity);
    double y = semi_major_axis * std::sqrt(1 - eccentricity * eccentricity) * std::sin(true_anomaly);

    vec3 p_hat = vec_eccentricity / eccentricity;
    vec3 q_hat = (h /h.norm()).prod(p_hat);

    x1 = p_hat * x + q_hat * y;
    x1 = x1 * m2 / (m1 + m2);
    x2 = p_hat * x + q_hat * y;
    x2 = x2 * -m1 / (m1 + m2);
}

void Euler::step() {
    compute_accelerations(a, x);

    for (int i = 0; i < m.size(); i++) {
        x[i] = x[i] + v[i] * dt;
        v[i] = v[i] + a[i] * dt / m[i];
    }
}

void EulerSwapped::step() {
    compute_accelerations(a, x);

    for (int i = 0; i < m.size(); i++) {
        v[i] = v[i] + a[i] * dt / m[i];
        x[i] = x[i] + v[i] * dt;
    }
}

void Leapfrog::step() {
    compute_accelerations(a, x);
    for (int i = 0; i < m.size(); i++) {
        x[i] = x[i] + v[i] * dt + a[i] * 0.5 * dt * dt;
    }

    std::vector<vec3> ap1(m.size());
    compute_accelerations(ap1, x);
    for (int i = 0; i < m.size(); i++) {
        v[i] = v[i] + (a[i] + ap1[i]) * 0.5 * dt / m[i];
    }
}

void RK2::step() {
    std::vector<vec3> a_star(m.size());
    compute_accelerations(a_star, x);
    for (int i = 0; i < m.size(); i++) {
        vec3 v_star = v[i] + a_star[i] * dt * theta / m[i];
        vec3 x_star = x[i] + v_star * dt * theta;
    }

    compute_accelerations(a, x);

    for (int i = 0; i < m.size(); i++) {
        x[i] = x[i] + v[i] * (1 - 1 / (2 * theta)) * dt + v[i] * 1 / (2 * theta) * dt;
        v[i] = v[i] + a[i] * dt / m[i];
    }
}

void RK4::step() {
    std::vector<vec3> k1(m.size());
    std::vector<vec3> k2(m.size());
    std::vector<vec3> k3(m.size());
    std::vector<vec3> k4(m.size());

    // tmp position vector
    std::vector<vec3> tmp(m.size());

    compute_accelerations(k1, x);
    for (int i = 0; i < m.size(); i++) {
        vec3 v1 = v[i] + k1[i] * 0.5 * dt / m[i];
        tmp[i] = x[i]  + v1 * dt;
    }

    compute_accelerations(k2, tmp);
    for (int i = 0; i < m.size(); i++) {
        vec3 v2 = v[i] + k2[i] * 0.5 * dt / m[i];
        tmp[i] = x[i]  + v2 * dt;
    }

    compute_accelerations(k3, tmp);
    for (int i = 0; i < m.size(); i++) {
        vec3 v3 = v[i] + k3[i] * dt / m[i];
        tmp[i] = x[i]  + v3 * dt;
    }

    compute_accelerations(k4, tmp);

    for (int i = 0; i < m.size(); i++) {
        x[i] = x[i] + v[i] * dt;
        v[i] = v[i] + (k1[i] + k2[i]*2 + k3[i]*2 + k4[i]) * dt / 6 / m[i];
    }
}

void RK4_md::step() {
    compute_acc_md(t, data, k1);
    for (int i = 0; i < data.extent(0); i++) {
        tmp[i] = data[i] + k1[i] * 0.5 * dt;
    }

    compute_acc_md(t + dt/2, tmp, k2);
    for (int i = 0; i < data.extent(0); i++) {
        tmp[i] = data[i] + k2[i] * 0.5 * dt;
    }

    compute_acc_md(t + dt/2, tmp, k3);
    for (int i = 0; i < data.extent(0); i++) {
        tmp[i] = data[i] + k3[i] * dt;
    }

    compute_acc_md(t + dt, tmp, k4);

    for (int i = 0; i < data.extent(0); i++) {
        data[i] = data[i] + (k1[i] + k2[i]*2 + k3[i]*2 + k4[i]) * dt / 6;
    }

    for (int i = 0; i < m.size(); i++) {
        x[i] = _data[i];
        v[i] = _data[m.size() + i];
    }
}


void RK45::step() {
    std::vector<vec3> k1(m.size());
    std::vector<vec3> k2(m.size());
    std::vector<vec3> k3(m.size());
    std::vector<vec3> k4(m.size());
    std::vector<vec3> k5(m.size());
    std::vector<vec3> k6(m.size());

    // tmp position vector
    std::vector<vec3> tmp(m.size());

    std::vector<vec3> e(m.size());

    double sum_e;

    do {
        compute_accelerations(k1, x);
        for (int i = 0; i < m.size(); i++) {
            vec3 v1 = v[i] + k1[i] * dt * B[1][0] / m[i];
            tmp[i] = x[i] + v1 * dt;
        }

        compute_accelerations(k2, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 v2 = v[i] + (k1[i] * B[2][0] + k2[i] * B[2][1]) * dt / m[i];
            tmp[i] = x[i] + v2 * dt;
        }

        compute_accelerations(k3, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 v3 = v[i] + (k1[i] * B[3][0] + k2[i] * B[3][1] + k3[i] * B[3][2]) * dt / m[i];
            tmp[i] = x[i] + v3 * dt;
        }

        compute_accelerations(k4, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 v4 = v[i] + (k1[i] * B[4][0] + k2[i] * B[4][1] + k3[i] * B[4][2] + k4[i] * B[4][3]) * dt / m[i];
            tmp[i] = x[i] + v4 * dt;
        }

        compute_accelerations(k5, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 v5 = v[i] + (k1[i] * B[5][0] + k2[i] * B[5][1] + k3[i] * B[5][2] + k4[i] * B[5][3] + k5[i] * B[5][4]) * dt / m[i];
            tmp[i] = x[i] + v5 * dt;
        }

        compute_accelerations(k6, tmp);

        sum_e = 0;
        for (int i = 0; i < m.size(); i++) {
            e[i] = k1[i] * (ch[0] - c[0]) + k2[i] * (ch[1] - c[1]) + k3[i] * (ch[2] - c[2]) + k4[i] * (ch[3] - c[3]) + k5[i] * (ch[4] - c[4]) + k6[i] * (ch[5] - c[5]);
            // physics ?
            sum_e += e[i].norm();
        }

        dt = 0.9 * dt * std::pow(eps / sum_e, 1.0/5.0);
    } while (sum_e > eps);

    for (int i = 0; i < m.size(); i++) {
        x[i] = x[i] + v[i] * dt;
        v[i] = v[i] + (k1[i] * ch[0] + k2[i] * ch[1] + k3[i] * ch[2] + k4[i] * ch[3] + k5[i] * ch[4] + k6[i] * ch[5]) * dt / m[i];
    }

    // TMP ??
    if (dt > base_dt) dt = base_dt * 1;
}

void DOP853::step() {
    using namespace dop853;

    std::vector<vec3>  k1(m.size());
    std::vector<vec3>  k2(m.size());
    std::vector<vec3>  k3(m.size());
    std::vector<vec3>  k4(m.size());
    std::vector<vec3>  k5(m.size());
    std::vector<vec3>  k6(m.size());
    std::vector<vec3>  k7(m.size());
    std::vector<vec3>  k8(m.size());
    std::vector<vec3>  k9(m.size());
    std::vector<vec3> k10(m.size());

    // tmp position vector
    std::vector<vec3> tmp(m.size());

    std::vector<vec3> e(m.size());

    // NOTE: ref impl takes an educated guess on the value of dt for step 0,
    // we're skipping that i'm lazy

    double err = 0;
    do {
        compute_accelerations(k1, x);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a21) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k2, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a31 + k2[i] * a32) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k3, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a41 + k3[i] * a43) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k4, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a51 + k3[i] * a53 + k4[i] * a54) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k5, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a61 + k4[i] * a64 + k5[i] * a65) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k6, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a71 + k4[i] * a74 + k5[i] * a75 + k6[i] * a76) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k7, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a81 + k4[i] * a84 + k5[i] * a85 + k6[i] * a86 + k7[i] * a87) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k8, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a91 + k4[i] * a94 + k5[i] * a95 + k6[i] * a96 + k7[i] * a97
                    + k8[i] * a98) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k9, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a101 + k4[i] * a104 + k5[i] * a105 + k6[i] * a106 
                    + k7[i] * a107 + k8[i] * a108 + k9[i] * a109) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        compute_accelerations(k10, tmp);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a111 + k4[i] * a114 + k5[i] * a115 + k6[i] * a116 
                    + k7[i] * a117 + k8[i] * a118 + k9[i] * a119 + k10[i] * a1110) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        // k2 reuse !! not a bug !! 
        compute_accelerations(k2, x);
        for (int i = 0; i < m.size(); i++) {
            vec3 _v = v[i] + (k1[i] * a121 + k4[i] * a124 + k5[i] * a125 + k6[i] * a126 
                    + k7[i] * a127 + k8[i] * a128 + k9[i] * a129 + k10[i] * a1210 + k2[i] * a1211) * dt / m[i];
            tmp[i] = x[i] + _v * dt;
        }

        // same for k3 !!1!
        compute_accelerations(k3, tmp);
        for (int i = 0; i < m.size(); i++) {
            k4[i] = k1[i] * b1 + k6[i] * b6 + k7[i] * b7 + k8[i] * b8 + k9[i] * b9 + k10[i] * b10 + k2[i] * b11 + k3[i] * b12;
            // might haved fucked up here (ref l. 651)
            // k5[i] = x[i] + (v[i] + k4[i] * dt / m[i]) * dt; 
            k5[i] = v[i] + k4[i] * dt / m[i]; 
        } 

        // error estimation
        // idk if all those vec3.norm() are good ideas
        err = 0;
        double err2 = 0;
        for (int i = 0; i < m.size(); i++) {
            double sk = a_tol + r_tol * std::max(std::abs(v[i].x), std::abs(k5[i].x));
            double err_i = (k4[i] - k1[i] * bhh1 - k9[i] * bhh2 - k3[i] * bhh3).x;
            err2 = err2 + (err_i / sk) * (err_i / sk);
            err_i = (k1[i] * er1 + k6[i] * er6 + k7[i] * er7 + k8[i] * er8 + k9[i] * er9 + k10[i] * er10 + k2[i] * er11 + k3[i] * er12).x;
            err = err + (err_i / sk) * (err_i / sk);

            sk = a_tol + r_tol * std::max(std::abs(v[i].y), std::abs(k5[i].y));
            err_i = (k4[i] - k1[i] * bhh1 - k9[i] * bhh2 - k3[i] * bhh3).y;
            err2 = err2 + (err_i / sk) * (err_i / sk);
            err_i = (k1[i] * er1 + k6[i] * er6 + k7[i] * er7 + k8[i] * er8 + k9[i] * er9 + k10[i] * er10 + k2[i] * er11 + k3[i] * er12).y;
            err = err + (err_i / sk) * (err_i / sk);

            sk = a_tol + r_tol * std::max(std::abs(v[i].z), std::abs(k5[i].z));
            err_i = (k4[i] - k1[i] * bhh1 - k9[i] * bhh2 - k3[i] * bhh3).z;
            err2 = err2 + (err_i / sk) * (err_i / sk);
            err_i = (k1[i] * er1 + k6[i] * er6 + k7[i] * er7 + k8[i] * er8 + k9[i] * er9 + k10[i] * er10 + k2[i] * er11 + k3[i] * er12).z;
            err = err + (err_i / sk) * (err_i / sk);
            std::cout << err << " " << err2 << std::endl;
        }
        double deno = err + 0.01 * err2;
        if (deno <= 0) deno = 1;
        // m.size() is supposed to be total number of dims (3 * m.size() ?)
        err = std::abs(dt) * err * std::sqrt(1.0 / (3 * m.size() * deno));

        // new dt
        // missing many parameters here (ref l.681)
        double fac11 = std::pow(err, 1.0 / 8.0 - beta * 0.2);
        double fac = fac11 / std::pow(facold, beta);
        std::cout << err << " " << fac << std::endl;
        fac = std::max(1.0/6.0, std::min(fac11, fac/0.9));
        dt = dt / fac;
        facold = std::max(err, 1e-4);
        std::cout << steps << " " << dt << " " << fac << " " << std::endl;
        // dt = 0.001;
        // err = 0.1;
    } while (err > 1.0);   

    // we skip both stiffness detection and dense output prep
    for (int i = 0; i < m.size(); i++) {
        x[i] = x[i] + v[i] * dt;
        v[i] = v[i] + k4[i] * dt / m[i];
    }
}

LinearMultistep::LinearMultistep(double target_t, int _back_steps, int bodies, int seed=0) : System{"LinearMultistep", target_t, bodies, seed}, back_steps{_back_steps} {
    prev.resize(m.size());

    for (int i = 0; i < m.size(); i++) {
        prev[i] = v[i];
    }

}

void LinearMultistep::step() {
    compute_accelerations(a, x);

    for (int i = 0; i < m.size() ;i++) {
        x[i] = x[i] + (v[i] * 1.5 - prev[i] * 0.5) * dt;
        prev[i] = v[i];
        v[i] = v[i] + a[i] * dt / m[i];
    }
}