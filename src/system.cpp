#include "system.hpp"

using std::chrono::high_resolution_clock;
using std::chrono::time_point;

config::config(int _bodies, fpoint_t _target_t, int _seed) {
    seed = _seed;
    bodies = _bodies;
    target_t = _target_t;

    m.resize(bodies);
    x.resize(bodies);
    v.resize(bodies);

    // generate random masses/positions/speeds
    std::srand(seed);
    for (int i = 0; i < bodies; i++) {
        m[i] = ((fpoint_t) (std::rand()) / RAND_MAX) * 2 + 3;
        x[i] = vec3{
            ((fpoint_t) (std::rand()) / RAND_MAX * 2 - 1) * 4,
            ((fpoint_t) (std::rand()) / RAND_MAX * 2 - 1) * 4,
            ((fpoint_t) (std::rand()) / RAND_MAX * 2 - 1) * 4,
        };
        v[i] = vec3{
            ((fpoint_t) (std::rand()) / RAND_MAX * 2 - 1),
            ((fpoint_t) (std::rand()) / RAND_MAX * 2 - 1),
            ((fpoint_t) (std::rand()) / RAND_MAX * 2 - 1),
        };
    }

    // offset those values so the barycenter is the origin of the coordinates
    vec3 barycenter_speed = vec3{0, 0, 0};
    vec3 barycenter_pos = vec3{0, 0, 0};
    fpoint_t m_sum = 0;
    for (int i = 0; i < bodies; i++) {
        barycenter_speed = barycenter_speed + v[i] * m[i];
        barycenter_pos = barycenter_pos + x[i] * m[i];
        m_sum += m[i];
    }
    barycenter_speed = barycenter_speed / m_sum;
    barycenter_pos = barycenter_pos / m_sum;

    for (int i = 0; i < bodies; i++) {
        v[i] = v[i] - barycenter_speed;
        x[i] = x[i] - barycenter_pos;
    }
}
