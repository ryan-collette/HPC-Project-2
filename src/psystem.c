#include "psystem.h"
#include "physics.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <stdio.h>

static Particle *particles; 
static size_t N_particles;

void ps_init(size_t p_count)
{
	free(particles);

	N_particles = p_count;
	particles = malloc(N_particles * sizeof(Particle));

	srand(time(0));
}

float randf()
{
	return (float)rand() / (float)RAND_MAX;
}

void ps_randomize(float radius)
{
	if (particles == NULL)
		return;
}

void ps_destroy()
{
	free(particles);
}
