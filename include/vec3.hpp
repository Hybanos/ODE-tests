#pragma once

#include <iostream>

#include "defs.hpp"

struct vec3 {
    fpoint_t x;
    fpoint_t y;
    fpoint_t z;

    inline vec3 operator-() {
        return vec3{-x, -y, -z};
    }

    inline vec3 operator+(vec3 v) {
        return vec3{x + v.x, y + v.y, z + v.z};
    }

    inline vec3 operator-(vec3 v) {
        return vec3{x - v.x, y - v.y, z - v.z};
    }

    inline vec3 operator*(vec3 v) {
        return vec3{x * v.x, y * v.y, z * v.z};
    }

    inline vec3 operator/(vec3 v) {
        return vec3{x / v.x, y / v.y, z / v.z};
    }

    inline vec3 operator+(fpoint_t d) {
        return vec3{x + d, y + d, z + d};
    }

    inline vec3 operator-(fpoint_t d) {
        return vec3{x - d, y - d, z - d};
    }

    inline vec3 operator*(fpoint_t d) {
        return vec3{x * d, y * d, z * d};
    }

    inline vec3 operator/(fpoint_t d) {
        return vec3{x / d, y / d, z / d};
    }

    inline fpoint_t dist_squared(vec3 o) {
        return (x - o.x) * (x - o.x) + (y - o.y) * (y - o.y) + (z - o.z) * (z - o.z);
    }

    inline fpoint_t dist(vec3 o) {
        return std::sqrt(dist_squared(o));
    }

    inline fpoint_t norm() {
        return std::sqrt(x * x + y * y + z * z);
    }

    inline vec3 unit() {
        return vec3{x, y, z} / norm();
    }

    inline vec3 prod(vec3 v) {
        return vec3{y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
    }

    inline fpoint_t dot(vec3 v) {
        return (*this * v).norm();
    }

    inline fpoint_t reduce() {
        return x + y + z;
    }

    void print() {
        std::cout << "vec3(" << x << "; " << y << "; " << z << ")" << std::endl;
    }
};