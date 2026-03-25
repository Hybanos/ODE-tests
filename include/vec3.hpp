#pragma once

#include <iostream>

struct vec3 {
    double x;
    double y;
    double z;

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

    inline vec3 operator+(double d) {
        return vec3{x + d, y + d, z + d};
    }

    inline vec3 operator-(double d) {
        return vec3{x - d, y - d, z - d};
    }

    inline vec3 operator*(double d) {
        return vec3{x * d, y * d, z * d};
    }

    inline vec3 operator/(double d) {
        return vec3{x / d, y / d, z / d};
    }

    inline double dist_squared(vec3 o) {
        return (x - o.x) * (x - o.x) + (y - o.y) * (y - o.y) + (z - o.z) * (z - o.z);
    }

    inline double dist(vec3 o) {
        return std::sqrt(dist_squared(o));
    }

    inline double norm() {
        return std::sqrt(x * x + y * y + z * z);
    }

    inline vec3 unit() {
        return vec3{x, y, z} / norm();
    }

    inline vec3 prod(vec3 v) {
        return vec3{y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
    }

    inline double dot(vec3 v) {
        return (*this * v).norm();
    }

    void print() {
        std::cout << "vec3(" << x << "; " << y << "; " << z << ")" << std::endl;
    }
};