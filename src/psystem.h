#pragma once
#include <stdlib.h>

extern double ps_G;
extern double ps_dt;
extern double ps_framerate;

size_t get_N_particles();

void ps_init(size_t p_count);
void ps_destroy();
void ps_run(double tspan);

void ps_randomize(double radius, double max_speed, double mass_min, double mass_max);

void ps_testcase();
