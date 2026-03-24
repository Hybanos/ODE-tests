#include "methods.hpp"

vec3 euler(std::function<vec3()>f, double dt) {
    vec3 dv = f();

    return dv * dt;
}