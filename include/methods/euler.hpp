#pragma once
#include <cmath>
#include <functional>

#include "vec3.hpp"

vec3 euler(std::function<vec3(void)>f, double dt);