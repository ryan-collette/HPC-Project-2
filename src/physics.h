#pragma once

#include <stdio.h>
#include <math.h>

typedef struct particle 
{
	double mass;
	double pos[3];
} Particle; 

const int N_PARTICLE_CMPTS = 4;

void gforce(double *pa, double *pb, double G, double *f)
{
	//minimum square distance between particles
	const double EPS = 1e-9;

	double dx = pb[0] - pa[0];
	double dy = pb[1] - pa[1];
	double dz = pb[2] - pa[2];
	double rr = dx*dx + dy*dy + dz*dz;
	double mag = G * pow(rr + EPS, -1.5);

	f[0] = dx * mag;	
	f[1] = dy * mag;
	f[2] = dz * mag;
}

/*
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
*/
