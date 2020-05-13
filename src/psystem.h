#pragma once
#include <stdlib.h>

void ps_init(size_t p_count);

void ps_destroy();

void ps_randomize(double radius, double max_speed, double mass_min, double mass_max);

void ps_run(double tspan);

