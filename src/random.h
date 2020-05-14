#include <stdlib.h>
#include <time.h>

const double PI = 3.14159265359;
const double PI_2 = 6.283185307;

void rand_init()
{
	srand(time(0));
}

double randf()
{
	return (double)rand() / (double)RAND_MAX;
}

double rand_frng(double min, double max)
{
	return min + (max - min) * randf();
}

void rand_sphere(double radius, double *p)
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
