#include "methods.hpp"

Exact::Exact(int max_steps) : System("Exact", max_steps) {
    // vec3 r_1 = (p2 - p1) * m2 / (m1 + m2);

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
    vec3 a = compute_acceleration(p1, p2);

    // speeds
    v1 = v1 + a * dt * m2;
    v2 = v2 - a * dt * m1; 

    // pos
    p1 = p1 + v1 * dt;
    p2 = p2 + v2 * dt;
}

void EulerSwapped::step() {
    vec3 a = compute_acceleration(p1, p2);

    // pos
    p1 = p1 + v1 * dt;
    p2 = p2 + v2 * dt;

    // speeds
    v1 = v1 + a * dt * m2;
    v2 = v2 - a * dt * m1; 
}

void Leapfrog::step() {
    vec3 a = compute_acceleration(p1, p2);

    p1 = p1 + v1 * dt + a * 0.5 * dt * dt * m2;
    p2 = p2 + v2 * dt - a * 0.5 * dt * dt * m1;

    vec3 ap1 = compute_acceleration(p1, p2);

    v1 = v1 + (a + ap1) * 0.5 * dt * m2;
    v2 = v2 - (a + ap1) * 0.5 * dt * m1;
}

void RK2::step() {
    // first step
    vec3 a_star = compute_acceleration(p1, p2);

    vec3 v1_star = v1 + a_star * dt * theta * m2;
    vec3 v2_star = v2 - a_star * dt * theta * m1; 

    vec3 p1_star = p1 + v1_star * dt * theta;
    vec3 p2_star = p2 + v2_star * dt * theta;

    // second step
    vec3 a = compute_acceleration(p1_star, p2_star);

    v1 = v1 + a * dt * m2;
    v2 = v2 - a * dt * m1; 

    p1 = p1 + v1_star * (1 - 1 / (2 * theta)) * dt + v1 * 1 / (2 * theta) * dt;
    p2 = p2 + v2_star * (1 - 1 / (2 * theta)) * dt + v2 * 1 / (2 * theta) * dt;
}

LinearMultistep::LinearMultistep(int max_steps, int _back_steps) : System{"LinearMultistep", max_steps}, back_steps{_back_steps} {

};

void LinearMultistep::step() {
    vec3 a = compute_acceleration(p1, p2);
    // speeds
    v1 = v1 + a * dt * m2;
    v2 = v2 - a * dt * m1; 

    // pos
    p1 = p1 + (v1 * 1.5 - prev1 * 0.5) * dt;
    p2 = p2 + (v2 * 1.5 - prev2 * 0.5) * dt;

    prev1 = v1;
    prev2 = v2;
}