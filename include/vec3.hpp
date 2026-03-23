struct vec3 {
    double x;
    double y;
    double z;

    // scalar
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

    // vector
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

    inline double dist_squared(vec3 v) {
        return (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z);
    }
};