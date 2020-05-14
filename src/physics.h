#pragma once

#include <stdio.h>
#include <math.h>

const double PI = 3.14159265359;
const double PI_2 = 6.283185307;

typedef struct particle 
{
	double m;
	double p[3];
	double v[3];
	double a[3];
} Particle; 

const Particle DEFAULT_PARTICLE = { 1, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }};
	 

void gforce(Particle *a, Particle *b, double G, double *f)
{
	//minimum square distance between particles
	const double EPS = 1e-9;

	double dx = b->p[0] - a->p[0];
	double dy = b->p[1] - a->p[1];
	double dz = b->p[2] - a->p[2];
	double rr = dx*dx + dy*dy + dz*dz;
	double mag = G * pow(rr + EPS, -1.5);

	f[0] = dx * mag;	
	f[1] = dy * mag;
	f[2] = dz * mag;
}

void write_particles(Particle *particles, int n, const char *fname)
{
	FILE *file = fopen(fname, "w");

	for (int i = 0; i < n; i++)
		fprintf(file, "%lf, %lf, %lf\n", particles[i].p[0], particles[i].p[1], particles[i].p[2]); 

	fclose(file);
}

void read_particles(Particle *particles, int n, const char *fname)
{
	FILE *file = fopen(fname, "r");

	for (int i = 0; i < n; i++)
		fscanf(file, "%lf %lf %lf", &particles[i].p[0], &particles[i].p[1], &particles[i].p[2]);

	fclose(file);
}
