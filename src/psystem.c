#include "psystem.h"
#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

double ps_G = 1.0;
double ps_dt = 0.001;
double ps_framerate = 10;

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

void ps_testcase()
{
	N_particles = 2;	

	Particle p1 = DEFAULT_PARTICLE;
	p1.p[0] = 1.0;
	p1.v[1] = 1.0;

	Particle p2 = DEFAULT_PARTICLE;
	p2.p[0] = -1.0;
	p2.v[1] = -1.0;

	particles[0] = p1;
	particles[1] = p2;	
}

static inline void ps_step()
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
				particles[i].a[k] += f[k] * particles[j].m;
				particles[j].a[k] -= f[k] * particles[i].m;
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
	char outfile[128]; 
	int n = (int)ceil(tspan / ps_dt);
	int frame_stride = ps_framerate > 0 ? (int)ceil(1.0 / (ps_dt * ps_framerate)) : -1;
	int prog_stride = n / 100;

	if (frame_stride > 0)
		write_particles(particles, N_particles, "frame0");

	for (int i = 1; i <= n; i++)
	{
		ps_step();

		if (i % prog_stride == 0)
		{
			double pct_prog = (double)i / n;
			printf("Progress: %.0f%%\n", pct_prog * 100);
		}

		if (frame_stride > 0 && i % frame_stride == 0)
		{
			snprintf(outfile, sizeof(outfile), "frame%d", i / frame_stride);
			write_particles(particles, N_particles, outfile);
		}
	}
}

