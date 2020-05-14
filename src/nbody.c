#include "psystem.h"
#include <stdio.h>

int main(int argc, char **argv)
{
	printf("Starting simulation.\n");

	size_t pcount = 100;	
	double tspan = 5.0;

	if (argc > 1)
		pcount = atoi(argv[1]);

	if (argc > 2)
		tspan = atof(argv[2]);

	if (argc > 4)
	{
		ps_G = atof(argv[3]);
		ps_dt = atof(argv[4]);
	}

	ps_init(pcount);

	printf("particle count = %ld\n", get_N_particles());
	printf("tspan = %f\n", tspan);
	printf("G = %f\n", ps_G);
	printf("dt = %f\n", ps_dt);

	ps_randomize(100, 1.0, 1.0, 1.0);
	ps_run(tspan);
	
	printf("Ending simulation.\n");
	ps_destroy();
	return 0;
}
