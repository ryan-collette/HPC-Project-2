#pragma once

#include <stdio.h>

const double PI = 3.14159265359;
const double PI_2 = 6.283185307;

typedef struct particle 
{
	double m;
	double x, y, z;
	double vx, vy, vz;
} Particle; 

void write_particles(Particle *arr, int n, const char *fname)
{
	FILE *file = fopen(fname, "w");

	for (int i = 0; i < n; i++)
		fprintf(file, "%lf, %lf, %lf\n", arr[i].x, arr[i].y, arr[i].z); 

	fclose(file);
}

void read_particles(Particle *arr, int n, const char *fname)
{
	FILE *file = fopen(fname, "r");

	for (int i = 0; i < n; i++)
	{
		fscanf(file, "%lf %lf %lf", &(arr+i)->x, &(arr+i)->y, &(arr+i)->z);
		printf("(%lf, %lf, %lf)\n", (arr+i)->x, (arr+i)->y, (arr+i)->z);
	}

	fclose(file);
}
