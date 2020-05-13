#include "psystem.h"
#include "physics.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

static Particle *particles; 
static size_t N_particles;

void ps_init(size_t p_count)
{
	free(particles);

	N_particles = p_count;
	particles = malloc(N_particles * sizeof(Particle));
	Particle new_part = { 1, 0, 0, 0, 0, 0, 0 };

	for (int i = 0; i < N_particles; i++)
		particles[i] = new_part; 

	srand(time(0));
}

double randf()
{
	return (double)rand() / (double)RAND_MAX;
}

double randf_rng(double min, double max)
{
	return min + (max - min) * randf();
}

void rnd_sphere(double radius, double *x, double *y, double *z)
{
	double r = radius * sqrt(randf()); 

	double theta = randf() * PI_2;
	double phi = randf() * PI_2;
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	double cos_phi = cos(phi);
	double sin_phi = sin(phi); 

	*x = r * cos_theta * sin_phi;
	*y = r * sin_theta * cos_phi; 
	*z = r * cos_phi; 
}

void ps_randomize(double radius, double max_speed, double mass_min, double mass_max)
{
	if (particles == NULL)
		return;

	double x, y, z;

	for (int i = 0; i < N_particles; i++)
	{
		(particles + i)->m = randf_rng(mass_min, mass_max);

 		rnd_sphere(radius, &x, &y, &z);
		(particles + i)->x = x;
		(particles + i)->y = y;
		(particles + i)->z = z;

		rnd_sphere(max_speed, &x, &y, &z);
		(particles + i)->vx = x;
		(particles + i)->vy = y;
		(particles + i)->vz = z;
	}

	write_particles(particles, N_particles, "test.csv");
	read_particles(particles, N_particles, "test.csv");
}

void ps_destroy()
{
	free(particles);
}
