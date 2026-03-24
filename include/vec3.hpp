struct vec3 {
    double x;
    double y;
    double z;

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
};
