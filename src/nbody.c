#include "psystem.h"
#include <stdio.h>
#include <sys/time.h>

static inline double wtime()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec / 1e6;
}

int main(int argc, char **argv)
{
	printf("Please enter any number of the following arguments:\n");
	printf("particle count, timespan, G (force constant), step size, framerate\n");

	size_t pcount = 2;
	double tspan = 5.0;

	scanf("%ld %lf %lf %lf %lf", &pcount, &tspan, &ps_G, &ps_dt, &ps_framerate);

	ps_init(pcount);
	ps_testcase();

	printf("\nStarting simulation with\n");
	printf("particle count = %ld\n", get_N_particles());
	printf("tspan = %lf\n", tspan);
	printf("G = %lf\n", ps_G);
	printf("dt = %lf\n", ps_dt);
	printf("framerate  = %lf\n", ps_framerate);

	//ps_randomize(100, 1.0, 1.0, 1.0);
	double t_start = wtime();
	ps_run(tspan);
	double t_end = wtime();
	
	printf("\nCompleted in %lf seconds\n", t_end - t_start);
	printf("\nEnding simulation.\n");
	ps_destroy();
	return 0;
}
