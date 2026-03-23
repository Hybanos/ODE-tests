struct vec3 {
    double x;
    double y;
    double z;

    inline double dist_squared(vec3 o) {
        return (x - o.x) * (x - o.x) + (y - o.y) * (y - o.y) + (z - o.z) * (z - o.z);
    }
};