#include "methods.hpp"

Exact::Exact(int max_steps, int bodies, int seed=0) : System("Exact", max_steps, bodies, seed) {

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
        v[i] = v[i] + (k1[i] + k2[i]*2 + k3[i]*3 + k4[i]) * dt / 6 / m[i];
    }
    
}

LinearMultistep::LinearMultistep(int max_steps, int _back_steps, int bodies, int seed=0) : System{"LinearMultistep", max_steps, bodies, seed}, back_steps{_back_steps} {
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