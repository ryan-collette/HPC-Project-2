#pragma once
#include <stdlib.h>
#include <stdbool.h>

extern double ps_G;
extern double ps_tspan;
extern double ps_dt;
extern double ps_framerate;

bool ps_is_scheduler();
size_t ps_N_particles();

void ps_init(size_t p_count);
void ps_destroy();
void ps_run();

void ps_randomize(double radius, double max_speed, double mass_min, double mass_max);
void ps_testcase();
