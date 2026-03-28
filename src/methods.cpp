#include "methods.hpp"

Exact::Exact(int max_steps, int bodies, int seed=0) : System("Exact", max_steps, bodies, seed) {

    if (m.size() > 2) {exit(1);}

    double &m1 = m[0];
    double &m2 = m[1];
    vec3 &p1 = x[0];
    vec3 &p2 = x[1];
    vec3 &v1 = v[0];
    vec3 &v2 = v[1];

    barycenter_pos = (p1 * m1 + p2 * m2) / (m1 + m2);
    barycenter_speed = (v1 * m1 + v2 * m2) / (m1 + m2);

    vec3 r = p1 - p2;
    vec3 v = v1 - v2;
    vec3 h = r.prod(v);
    double mu = m1 * m2 / (m1 + m2);

    double eps = (v.norm() * v.norm()) / 2 - gamma * (m1 + m2) / r.norm();

    semi_major_axis = - gamma * (m1 + m2) / 2 * eps;
    vec3 vec_eccentricity = v.prod(h) / gamma * (m1 + m2) - r / r.norm();
    eccentricity = vec_eccentricity.norm();

    double cos_nu = vec_eccentricity.dot(r) / (eccentricity * r.norm());
    double nu = std::acos(cos_nu);

    if (r.dot(v) < 0) nu = 2*3.141592 - nu;

    double tan_E_0_2 = std::sqrt((1 - eccentricity) / (1 + eccentricity)) * std::tan(nu / 2);
    E_0 = std::atan(tan_E_0_2) * 2;

    mean_movement = std::sqrt(gamma * (m1 + m2) / semi_major_axis * semi_major_axis * semi_major_axis);

    std::cout << semi_major_axis << std::endl;
}

vec3 Exact::compute_pos(double t) {
    double M = mean_movement * t;
    
    return vec3{};
}

vec3 Exact::compute_barycenter(double t) {
    return barycenter_pos + barycenter_speed * t;
}

void Exact::step() {

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