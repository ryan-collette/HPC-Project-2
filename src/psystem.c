#include "psystem.h"
#include "physics.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

double ps_G = 1.0;
double ps_dt = 0.0001;
int ps_framestride = 10;

static Particle *particles; 
static size_t N_particles;

size_t get_N_particles()
{
	return N_particles;
}

void ps_init(size_t p_count)
{
	free(particles);

	N_particles = p_count;
	particles = malloc(N_particles * sizeof(Particle));

	for (int i = 0; i < N_particles; i++)
		particles[i] = DEFAULT_PARTICLE; 

	srand(time(0));
}

void ps_destroy()
{
	free(particles);
}

double randf()
{
	return (double)rand() / (double)RAND_MAX;
}

double randf_rng(double min, double max)
{
	return min + (max - min) * randf();
}

void rnd_sphere(double radius, double *p)
{
	double r = radius * sqrt(randf()); 

	double theta = randf() * PI_2;
	double phi = randf() * PI_2;
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	double cos_phi = cos(phi);
	double sin_phi = sin(phi); 

	p[0] = r * cos_theta * sin_phi;
	p[1] = r * sin_theta * cos_phi; 
	p[2] = r * cos_phi; 
}

void ps_randomize(double radius, double max_speed, double mass_min, double mass_max)
{
	if (particles == NULL)
		return;

	for (int i = 0; i < N_particles; i++)
	{
		(particles + i)->m = randf_rng(mass_min, mass_max);
 		rnd_sphere(radius, particles[i].p);
		rnd_sphere(max_speed, particles[i].v);
	}
}

static void ps_step()
{
	double h = ps_dt * 0.5;
	double f[3];

	for (int i = 0; i < N_particles; i++)
		for (int k = 0; k < 3; k++)
			particles[i].a[k] = 0.0;

	for (int i = 0; i < N_particles - 1; i++)
		for (int j = i+1; j < N_particles; j++)
		{
			gforce(particles + i, particles + j, ps_G, f);		
			for (int k = 0; k < 3; k++)
			{
				particles[i].a[k] += f[k];
				particles[j].a[k] -= f[k];
			}	
		}

	for (int i = 0; i < N_particles; i++)
		for (int k = 0; k < 3; k++)
		{
			double v = particles[i].v[k];
			double a = particles[i].a[k];
			double vhalf = v + h * a;
			particles[i].p[k] += h * vhalf;
			particles[i].v[k] = vhalf + h * a;
		}
}

void ps_run(double tspan)
{
	int n = (int)ceil(tspan / ps_dt);
	char outfile[64]; 

	write_particles(particles, N_particles, "frame0");

	for (int i = 1; i <= n; i++)
	{
		ps_step();

		if (i % ps_framestride == 0)
		{
			snprintf(outfile, sizeof(outfile), "frame%d", i / ps_framestride);
			write_particles(particles, N_particles, outfile);
		}
	}
}

