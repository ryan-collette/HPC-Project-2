#include "psystem.h"
#include <stdio.h>

int main(int argc, char **argv)
{
	printf("Starting simulation.\n");
	ps_init(10000);

	ps_randomize(1.0);
	
	printf("Ending simulation.\n");
	ps_destroy();
	return 0;
}
