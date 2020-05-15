#include "psystem.h"
#include "physics.h"
#include "random.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>

double ps_G = 1.0;
double ps_tspan = 1.0;
double ps_dt = 0.001;
double ps_framerate = 10;

typedef struct segment
{
	int index, start, count;
	int owner_rank;
	bool current, held;
} Segment;

typedef struct block
{
	int owner_rank;
	bool current, held;
} Block;

typedef enum jobType { NO_JOB, STEP, UPDATE_FORCES } JobType;

static const int SCHEDULER_RANK = 0;

static int rank, N_procs;

static Segment *segments;
static int N_segments;

static Block *blocks;
static int N_blocks;

static Particle *particles; 
static size_t N_particles;

static void init_scheduler()
{
	N_segments = N_procs;
	segments = malloc(N_segments * sizeof(Segment));

	Segment s = { 0, 0, 0, rank, true, false }; 
	int r = N_particles % N_segments;
	int count = N_particles / N_segments;
	for (int i = 0; i < N_segments; i++)
	{
		s.index = i;

		if (i < r)
		{
			s.count = count + 1;
			s.start = s.count * i;
		}
		else
			s.count = count;
			s.start = (count + 1) * r + count * (i - r);

		segments[i] = s;
	}

	N_blocks = N_segments * N_segments;
	blocks = malloc(N_blocks * sizeof(Block));

	Block b = { rank, false, false }; 
	
	for (int r = 0; r < N_segments; r++)
		for (int c = 0; c < N_segments; c++)
			blocks[r * N_segments + c] = b;

	printf("N_procs: %d\n", N_procs);

	printf("Segments:\n");
	for (int i = 0; i < N_segments; i++)
		printf("{ %d, %d, %d }\n", segments[i].index, segments[i].start, segments[i].count);
}

static JobType request_job(Segment *row, Segment *col)
{
	for (int r = 0; r < N_segments; r++)
	{
		int c = -1;
		for (int i = r; i < N_segments; i++)
			if (!blocks[r * N_segments + i].current)
			{
				c = i;
				break;
			}

		if (c < 0)
		{
			if (segments[r].current)
			{
				for (int i = 0; i < N_segments; i++)
				{
					blocks[r * N_segments + i].current = false;
					blocks[i * N_segments + r].current = false; 
				} 

				return NO_JOB;
			}
			else
			{
				return STEP;
			}	
		}
		else
		{
			return UPDATE_FORCES;
		}
	}	

	return NO_JOB;
}

bool ps_is_scheduler()
{
	return rank == SCHEDULER_RANK;
}

size_t ps_N_particles()
{
	return N_particles;
}

void ps_init(size_t p_count)
{
	N_particles = p_count;
	particles = malloc(N_particles * sizeof(Particle));

	for (int i = 0; i < N_particles; i++)
		particles[i] = DEFAULT_PARTICLE; 

	rand_init();

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &N_procs);

	if (rank == SCHEDULER_RANK)
		init_scheduler();
}

void ps_destroy()
{
	MPI_Finalize();
	free(particles);
	free(segments);
	free(blocks);
}

void ps_randomize(double radius, double max_speed, double mass_min, double mass_max)
{
	for (int i = 0; i < N_particles; i++)
	{
		(particles + i)->m = rand_frng(mass_min, mass_max);
 		rand_sphere(radius, particles[i].p);
		rand_sphere(max_speed, particles[i].v);
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

void ps_run()
{
	char outfile[128]; 
	int n = (int)ceil(ps_tspan / ps_dt);
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

