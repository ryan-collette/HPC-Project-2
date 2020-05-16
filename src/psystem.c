#include "psystem.h"
#include "physics.h"
#include "random.h"

#include <stdbool.h>
#include <mpi.h>

double ps_G = 1.0;
double ps_tspan = 1.0;
double ps_dt = 0.001;
double ps_framerate = 0;

typedef struct job
{
	//two regions of particles for force-pair calcultation
	int row, col;

	MPI_Request request;
} Job;

static const int SCHEDULER_RANK = 0;
static const int BLOCKS_PER_PROC = 1;
static const int ROW_TAG = 200;
static const int COL_TAG = 201;
static const int ACCEL_TAG = 300; 

static int rank, N_procs;

//domain is divided into N_segs pieces called segments
static int N_segs;

static int parts_per_seg;  
static int cmpts_per_seg;

/********* scheduler data *********/

static Particle *particles;
static int N_particles;

static double *vels;
static int N_vels;

static double *accels; 
static int N_accels; 

//calculations from worker threads are stored here
//before being accumulated in scheduler data 
static double *work_buffer;
static int N_work_buffer; 

/*********************************/

/********* worker data  *********/

//workers perform force pair calculations on these
static Particle *particle_block;
static int N_particle_block;

static double *accel_block;
static int N_accel_block;

/********************************/

//flag used to stop worker procs
static bool terminate = false; 

static void run_scheduler();
static void run_worker(); 

size_t ps_N_particles()
{
	return N_particles;
}

bool ps_is_scheduler()
{
	return rank == SCHEDULER_RANK;
}

void ps_init(size_t p_count)
{
	rand_init();

	MPI_Init(NULL, NULL);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &N_procs);

	N_segs = N_procs * BLOCKS_PER_PROC; 

	p_count += p_count % N_segs;	
	parts_per_seg = p_count / N_segs;
	cmpts_per_seg = parts_per_seg * 3;

	if (ps_is_scheduler())
	{
		N_particles = p_count;
		N_vels = N_particles * 3;
		N_accels = N_particles * 3;
		N_work_buffer = p_count * 6; 

		particles = malloc(N_particles * sizeof(Particle));
		vels = malloc(N_vels * sizeof(double));
		accels = malloc(N_accels * sizeof(double));
		work_buffer = malloc(N_work_buffer * sizeof(double));	

		Particle p = { 1, { 0, 0, 0}};
		for (int i = 0; i < N_particles; i++)
			particles[i] = p;
	}	
	else
	{
		N_particle_block = parts_per_seg * 2;
		N_accel_block = cmpts_per_seg * 2;
	
		particle_block = malloc(N_particle_block * sizeof(Particle));
		accel_block = malloc(N_accel_block * sizeof(double));	
	}
}

void ps_destroy()
{
	free(particles);
	free(vels);
	free(accels);
	free(work_buffer);
	free(particle_block);
	free(accel_block);

	MPI_Finalize();
}

void ps_run()
{
	if (ps_is_scheduler())
		run_scheduler();
	else
		run_worker();
}

void ps_randomize(double radius, double max_speed,
                  double mass_min, double mass_max)
{
	for (int i = 0; i < N_particles; i++)
	{
		particles[i].mass = rand_frng(mass_min, mass_max);
		rand_sphere(radius, particles[i].pos);
		rand_sphere(max_speed, vels + 3*i);
	}
}

void ps_testcase()
{
	Particle p1 = { 1, { 1, 0, 0 } };
	vels[0] = 1;
	
	Particle p2 = { 1, { -1, 0, 0 } };
	vels[4] = -1;

	particles[0] = p1;
	particles[1] = p2;
}

//adds accelerations from worker threads to accel array
static void merge_jobs(Job *pending, bool *current)
{
	for (int i = 1; i < N_procs; i++)
	{
		int flag;
		Job job = pending[i];
		MPI_Test(&job.request, &flag, MPI_STATUS_IGNORE);

		if (flag && job.row > 0 && job.col > 0)
		{
			double *row_src = work_buffer + i * N_accel_block;
			double *col_src = row_src + cmpts_per_seg;
			double *row_dest = accels + job.row * cmpts_per_seg;
			double *col_dest = accels + job.col * cmpts_per_seg;

			for (int j = 0; j < cmpts_per_seg; j++)
			{
				row_dest[j] += row_src[j];
				col_dest[j] += col_src[j];
			} 

			current[job.row * N_segs + job.col] = true;
			current[job.col * N_segs + job.row] = true;

			pending[i].row = -1;
			pending[i].col = -1;
		}
	}
}

//Attempts to queue a new job for the given block.
//If successful, merges the completed job from the chosen worker proc.
static bool try_queue_job(int r, int c, Job *pending, bool *current)
{
	int rank = -1;

	//find free worker proc
	for (int i = 1; i < N_procs; i++)
		if (pending[i].row < 0 || pending[i].col < 0)
		{
			rank = i;
			break;
		}

	if (rank < 0)
		return false;

	Job job = { .row = r, .col = c };

	Particle *row = particles + r * parts_per_seg;
	MPI_Isend(row, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
	          rank, ROW_TAG, MPI_COMM_WORLD, &job.request);	

	Particle *col = particles + c * parts_per_seg;
	MPI_Isend(col, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
	          rank, COL_TAG, MPI_COMM_WORLD, &job.request);

	MPI_Irecv(work_buffer + rank * N_accel_block, N_accel_block, MPI_DOUBLE,
	          rank, ACCEL_TAG, MPI_COMM_WORLD, &job.request); 
	
	pending[rank] = job;
	return true;
} 

//scan through non-current blocks and queue jobs
static void update_jobs(Job *pending, bool *current)
{
	static int i = 0;

	merge_jobs(pending, current);

	if (i == N_segs * N_segs)
	{
		if (current[0])
			i = 0;
		else
			return;
	}

	while (i < N_segs)
	{
		int r = i / N_segs;
		int c = i % N_segs;

		//no workers available
		if (!current[i] && !try_queue_job(r, c, pending, current))
			return;

		i++;
	}
}

//perform integration step for segment i
static void step(int index)
{
	double h = 0.5 * ps_dt;

	Particle *seg = particles + index * parts_per_seg;
	double *seg_vels = vels + index * cmpts_per_seg;	
	double *seg_accels = accels + index * cmpts_per_seg;	

	for (int i = 0 ; i < parts_per_seg; i++)
	for (int k = 0; k < 3; k++)
	{
		double v = seg_vels[i*3 + k];
		double a = seg_accels[i*3 + k];
		double vhalf = v + h * a;
		seg[i].pos[k] += h * vhalf;
		seg_vels[i*3 + k] = vhalf + h * a;	
	}	
}

static void run_scheduler()
{
	int n = (int)ceil(ps_tspan / ps_dt);
	int frame_stride = ps_framerate > 0 ? (int)ceil(1.0 / (ps_dt * ps_framerate)) : -1;
	int prog_stride = n / 100;

	bool *current = malloc(N_segs * N_segs * sizeof(bool));
	Job *pending = malloc(N_procs * sizeof(Job)); 
	int seg_index = 0;

	while (n > 0)
	{
		update_jobs(pending, current);

		bool row_current = true;
		for (int i = 0; i < N_segs && row_current; i++)
			row_current &= current[seg_index * N_segs + i]; 

		if (row_current)
		{
			step(seg_index);

			for (int i = 0; i < N_segs; i++)
				current[seg_index * N_segs + i] = false;

			seg_index++;
			if (seg_index == N_procs)
			{
				seg_index = 0;
				n--;
			}	
		}
	}

	free(current);
	free(pending);
}

static void tick_worker()
{
	Particle *row = particle_block;
	Particle *col = particle_block + parts_per_seg;
	double *row_accels = accel_block;
	double *col_accels = accel_block + cmpts_per_seg;
	double f[3];

	MPI_Recv(row, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
	         SCHEDULER_RANK, ROW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
	MPI_Recv(col, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
	         SCHEDULER_RANK, COL_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	for (int i = 0; i < N_accel_block; i++)
		accel_block[i] = 0.0;		
	
	for (int i = 0; i < parts_per_seg; i++)
	{
		double *r_accel = row_accels + 3*i;
		
		for (int j = 0; j < parts_per_seg; j++)
		{
			double *c_accel = col_accels + 3*j;
			gforce(row[i].pos, col[j].pos, ps_G, f);

			for (int k = 0; k < 3; k++)
			{
				r_accel[k] += f[k] * col[j].mass;
				c_accel[k] -= f[k] * row[i].mass;
			}
		}
	} 

	MPI_Send(accel_block, N_accel_block, MPI_DOUBLE,
	         SCHEDULER_RANK, ACCEL_TAG, MPI_COMM_WORLD);
}

static void run_worker()
{
	while (!terminate)
		tick_worker();
} 

