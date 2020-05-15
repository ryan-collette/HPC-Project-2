#include "psystem.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

static inline double wtime()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec / 1e6;
}

static void print_help()
{
	printf("Please execute the program with any number of the following command line arguments in order:\n");
	printf("-particle count\n");
	printf("-timespan: how long the simulation runs\n");
	printf("-G: gravitational constant\n");
	printf("-dt: integration step size\n");
	printf("-framerate:\n");
	printf("--outputs each frame as a .csv file at the given rate\n");
	printf("--outputs nothing for framerate <= 0.0\n");
}

static void print_params()
{
	printf("\nStarting simulation with:\n");
	printf("particle count = %ld\n", ps_N_particles());
	printf("tspan = %lf\n", ps_tspan);
	printf("G = %lf\n", ps_G);
	printf("dt = %lf\n", ps_dt);
	printf("framerate  = %lf\n\n", ps_framerate);
}

int main(int argc, char **argv)
{
	if (argc > 1 && strcmp(argv[1], "help") == 0)
	{
		print_help();
		return 0;
	}

	size_t pcount = 2;
	if (argc > 1)
	{
		pcount = atoi(argv[1]);
		pcount = pcount < 1 ? 1 : pcount;
	}

	ps_init(pcount);

	if (argc > 2)
		ps_tspan = atof(argv[2]);
	if (argc > 3)
		ps_G = atof(argv[3]);
	if (argc > 4)
		ps_dt = atof(argv[4]);
	if (argc > 5)
		ps_framerate = atof(argv[5]);	

	if (ps_is_scheduler())
	{
		ps_randomize(100, 1.0, 1.0, 1.0);

		print_params();

		double t_start = wtime();
		//ps_run(tspan);
		double t_end = wtime();

		printf("\nCompleted in %lf seconds\n", t_end - t_start);
		printf("\nEnding simulation.\n");
	}
	else
		;
		//ps_run(tspan);

	ps_destroy();
	return 0;
}
