#include "psystem.h"
#include "physics.h"
#include "random.h"

#include <stdbool.h>
#include <mpi.h>
#include <assert.h>
#include <time.h>

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

	assert(N_procs > 1);

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

//Adds accelerations from worker threads to accel array
//Returns number of jobs still pending
static int merge_jobs(Job *jobs, bool *current)
{
	int count = 0;

	for (int i = 1; i < N_procs; i++)
	{
		if (jobs[i].row < 0 || jobs[i].col < 0) //job not pending
			continue;

		int complete = false;
		MPI_Test(&jobs[i].request, &complete, MPI_STATUS_IGNORE);
		int r_ind = jobs[i].row * N_segs + jobs[i].col;
		int c_ind = jobs[i].col * N_segs + jobs[i].row;

		if (complete && !current[r_ind] && !current[c_ind])
		{
			printf("Merging job (%d, %d) from %d\n", jobs[i].row, jobs[i].col, i);
/*
			double *row_src = work_buffer + i * N_accel_block;
			double *col_src = row_src + cmpts_per_seg;
			double *row_dest = accels + jobs[i].row * cmpts_per_seg;
			double *col_dest = accels + jobs[i].col * cmpts_per_seg;

			for (int j = 0; j < cmpts_per_seg; j++)
			{
				row_dest[j] += row_src[j];
				col_dest[j] += col_src[j];
			} 

*/
			current[r_ind] = true;
			current[c_ind] = true;

			jobs[i].row = -1;
			jobs[i].col = -1;
		}
		else
			count++;
	}

	return count;
}

//Queues the Job corresonding to job_i if a free worker is available.
//Returns true if successful, flase otherwise.
static bool try_queue_job(int job_i, Job *jobs, bool *current)
{
	int r = job_i / N_segs;
	int c = job_i % N_segs;

	for (int i = 1; i < N_procs; i++)
		if (jobs[i].row == r && jobs[i].col == c ||
		    jobs[i].col == r && jobs[i].row == c)
			return true;

	int rank = -1;

	//find free worker
	for (int i = 1; i < N_procs; i++)
		if (jobs[i].row < 0 || jobs[i].col < 0)
		{
			rank = i;
			break;
		}

	if (rank < 0)
		return false;

	Particle *row = particles + r * parts_per_seg;
	Particle *col = particles + c * parts_per_seg;
	MPI_Request row_request, col_request, accel_request;

	MPI_Isend(row, 0, MPI_DOUBLE, rank, ROW_TAG, MPI_COMM_WORLD, &row_request);	

	MPI_Isend(col, 0, MPI_DOUBLE, rank, COL_TAG, MPI_COMM_WORLD, &col_request);

	MPI_Irecv(work_buffer, 0, MPI_DOUBLE, rank, ACCEL_TAG, MPI_COMM_WORLD, &accel_request); 

	//MPI_Isend(row, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
	//          rank, ROW_TAG, MPI_COMM_WORLD, &row_request);	

	//Particle *col = particles + c * parts_per_seg;
	//MPI_Isend(col, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
	//          rank, COL_TAG, MPI_COMM_WORLD, &col_request);

	//MPI_Irecv(work_buffer + rank * N_accel_block, N_accel_block, MPI_DOUBLE,
	//          rank, ACCEL_TAG, MPI_COMM_WORLD, &accel_request); 

	jobs[rank].row = r;
	jobs[rank].col = c;
	jobs[rank].request = accel_request;
	printf("queuing job: (%d, %d) on %d\n", r, c, rank);
	return true;
} 

//perform integration step for segment i
static bool try_step(int seg_i, bool *current)
{
	bool row_current = true;
	for (int i = 0; i < N_segs && row_current; i++)
		row_current &= current[seg_i * N_segs + i]; 

	if (!row_current)
		return false;
/*
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
*/

	//row of blocks can now be computed for the next step
	for (int i = 0; i < N_segs; i++)
		current[seg_i * N_segs + i] = false;

	return true;
}

static void run_scheduler()
{
	Job *jobs = malloc(N_procs * sizeof(Job)); 
	Job null_job = { .row = -1, .col = -1, .request = MPI_REQUEST_NULL };
	for (int i = 0; i < N_procs; i++)
		jobs[i] = null_job;

	int N_current = N_segs * N_segs;
	//which blocks have completed for the current step?
	bool *current = malloc(N_current * sizeof(bool));
	for (int i = 0; i < N_current; i++)
		current[i] = false;

	int step_count = (int)ceil(ps_tspan / ps_dt);
	int frame_stride = ps_framerate > 0 ?
		(int)ceil(1.0 / (ps_dt * ps_framerate)) : -1;

	int seg_i = 0;
	int job_i = 0;
	int job_count = N_current * step_count;	

	while(step_count > 0)
	{
		merge_jobs(jobs, current);

		if (try_step(seg_i, current))
		{
			seg_i++;
			if (seg_i == N_segs)
			{
				seg_i = 0;
				step_count--;
			}	
		}

		while (job_count > 0 &&
			  (current[job_i] || try_queue_job(job_i, jobs, current)))
		{
			job_i = (job_i + 1) % N_current;
			job_count--;
		}

		double x = 0;
		for (int i = 0; i < 1000000; i++)
			x *= i;

		for (int r = 0; r < N_segs; r++)
		{
			for (int c = 0; c < N_segs; c++)
				printf("%d ", current[r * N_segs + c]);	
			printf("\n");
		}	
	}

	int dummy[0];
	MPI_Request stop_request;
	MPI_Ibcast(dummy, 0, MPI_INT, SCHEDULER_RANK, MPI_COMM_WORLD, &stop_request); 

	free(current);
	free(jobs);
}

static void compute_block(Particle *row, Particle *col,
                          double *row_accels, double *col_accels)
{
	/*
	double f[3];
	
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
	*/
}

static void tick_worker(MPI_Request *stop_request, int *terminate)
{
	Particle *row = particle_block;
	Particle *col = particle_block + parts_per_seg;
	double *row_accels = accel_block;
	double *col_accels = accel_block + cmpts_per_seg;

	MPI_Request row_request, col_request;
	int row_ready = false, col_ready = false;

	MPI_Irecv(row, 0, MPI_DOUBLE, SCHEDULER_RANK, ROW_TAG, MPI_COMM_WORLD, &row_request);	
	MPI_Irecv(col, 0, MPI_DOUBLE, SCHEDULER_RANK, COL_TAG, MPI_COMM_WORLD, &col_request);

//	MPI_Irecv(row, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
//	          SCHEDULER_RANK, ROW_TAG, MPI_COMM_WORLD, &row_request);	
//	MPI_Irecv(col, parts_per_seg * N_PARTICLE_CMPTS, MPI_DOUBLE,
 //	          SCHEDULER_RANK, COL_TAG, MPI_COMM_WORLD, &col_request);

	while(!row_ready || !col_ready)
	{
		MPI_Test(&row_request, &row_ready, MPI_STATUS_IGNORE);
		MPI_Test(&col_request, &col_ready, MPI_STATUS_IGNORE);
		MPI_Test(stop_request, terminate, MPI_STATUS_IGNORE);

		if (*terminate)
			return;
	}

	compute_block(row, col, row_accels, col_accels);

	MPI_Send(accel_block, 0, MPI_DOUBLE, SCHEDULER_RANK, ACCEL_TAG, MPI_COMM_WORLD);

	//MPI_Send(accel_block,_accel_block, MPI_DOUBLE,
	//         SCHEDULER_RANK, ACCEL_TAG, MPI_COMM_WORLD);
}

static void run_worker()
{
	int dummy[0];
	MPI_Request stop_request;
	int terminate = false;

	MPI_Ibcast(dummy, 0, MPI_INT, SCHEDULER_RANK, MPI_COMM_WORLD, &stop_request); 

	while (!terminate)
	{
		tick_worker(&stop_request, &terminate);
		MPI_Test(&stop_request, &terminate, MPI_STATUS_IGNORE);
	}
} 
