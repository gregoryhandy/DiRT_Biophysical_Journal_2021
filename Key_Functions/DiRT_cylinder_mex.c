/*
** The corresponding mex code for the DiRT simultion in the cylinder
** All parameters are defined in the correspond Matlab files
**
** Model details found in
** Handy G, Lawley SD, Revising Berg-Purcell for finite receptor kinetics,
** Biophysical Journal (2021), doi: https://doi.org/10.1016/j.bpj.2021.03.021.
** 
** Written by Gregory Handy, 04/05/2021
*/
#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"


/* returns a pseudorandom value drawn from the standard normal distribution */
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


// Receptor/Capture Region Structure
struct Receptor
{
	double center_x;
	double center_y; 
	double radius;
	//double start_pos;  // start of receptor
	//double end_pos;    // end of receptor
	double switch_time; // time until it switches back to the absorbing state
	double status; // 0 if open, 1 if occupied
};

// Particle Structure
struct Particle
{
	//current position
	double x_pos2; // x position
	double y_pos2; // y position
	double z_pos2; // z position
	
	// previous position
	double x_pos1; // x position
	double y_pos1; // y position
	double z_pos1; // z position
	
	int status;   // 0 for diffusing in cleft, -1 for diffusing left of the cleft, 1 for diffusing right of the cleft, 2 for captured/escaped
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double max_time, dt;
	double D, R, L_z;
	double x_start_loc, y_start_loc, z_start_loc; 
	double tau_r;
	int N, num_trials, M, print_num, max_time_points;
	
	double *rec_centers, rec_rad, start_std;
				
   /******
    * Import variables from matlab
    * This is messy looking and is specific to mex.
    * Ignore if you're implementing this outside of mex.
    *******/
	max_time =  mxGetScalar(prhs[0]);
	dt = mxGetScalar(prhs[1]);
	D =  mxGetScalar(prhs[2]);
	R =  mxGetScalar(prhs[3]);
	L_z =  mxGetScalar(prhs[4]);
	x_start_loc =  mxGetScalar(prhs[5]);
	y_start_loc =  mxGetScalar(prhs[6]);
	z_start_loc =  mxGetScalar(prhs[7]);
	N = (int)mxGetScalar(prhs[8]);
	tau_r =  mxGetScalar(prhs[9]);
	M = (int)mxGetScalar(prhs[10]);
	
	int m1, m2;
	
	int total_prints = 0;
	
	rec_centers = mxGetPr(prhs[11]);
	m1 = mxGetM(prhs[11]);
	m2 = mxGetN(prhs[11]);
	if(m1!=M && m2 != 2){
	    mexErrMsgTxt("rec_centers should be M x 2");
	}
	
	rec_rad = mxGetScalar(prhs[12]);
	start_std = mxGetScalar(prhs[13]);
	print_num = (int) mxGetScalar(prhs[14]);
	max_time_points = (int) mxGetScalar(prhs[15]); 

	struct Receptor* receptors = mxMalloc(M*sizeof(struct Receptor));
	
	// allocate memory for particles
	struct Particle* particles = mxMalloc(N*sizeof(struct Particle));
	
	
	double k_constant = sqrt(2*D*dt);
	
	// trial initialization
	double current_time = 0;
	// variables of interest
	// particles_remaining_v1: particle is removed from the domain once the receptor recharges
	int particles_remaining_v1 = N; 
	// particles_remaining_v2: particle is removed from the domain once it is captured by a receptor
	int particles_remaining_v2 = N;
	int total_captured_v1 = 0;
	int total_captured_v2 = 0;
	int num_recs_available = M;
	// helps keep track of when to store the results 
	int time_count = 1;
		
	// loop variables (could be indeclared at the start of each loop)
	int j,h,rec_loop,outer_count;
	h = 0;
		
	// used for particle reflection
	double x1, x2, y1, y2;
	double m,a,b,c;
	double x_inter1, y_inter1, x_inter2, y_inter2;
	double d1, d2, d;
	double x_inter_final, y_inter_final;
	double tangent_slope, C;
	double x_reflected, y_reflected;
	
	
	/* Allocate output vector */
	double *output_matrix_v1;
	plhs[0] = mxCreateDoubleMatrix(max_time_points,4, mxREAL);
	output_matrix_v1 = mxGetPr(plhs[0]);
	
	double *output_matrix_v2;
	plhs[1] = mxCreateDoubleMatrix(max_time_points,4, mxREAL);
	output_matrix_v2 = mxGetPr(plhs[1]);
	
	
	// all receptors start in the available state
	for(j=0;j<M;j++){
		receptors[j].radius = rec_rad;
		
		receptors[j].center_x = rec_centers[j];
		receptors[j].center_y = rec_centers[j + M];
		
		receptors[j].switch_time = -10;
		receptors[j].status = 0;
	}
	
	// seed random number generator
	srand(time(NULL));
	
	// particles starting position
	for(j=0;j<N;j++){
		particles[j].x_pos1 = randn(x_start_loc,start_std);
		particles[j].y_pos1 = randn(y_start_loc,start_std);
		particles[j].z_pos1 = randn(z_start_loc,start_std);
		particles[j].status = 0;
	}
			
	// main loop
	while((particles_remaining_v2 > 0 || num_recs_available < M) && current_time <= max_time){

		// Update the status of the receptors
		for(rec_loop=0;rec_loop<M;rec_loop++)
		{
			// receptor switchs from closed (1) to open (0)
			// Particle is officially counted as "captured"
	        if(receptors[rec_loop].status == 1 && current_time>receptors[rec_loop].switch_time)
			{
				// particle has been captured and is removed from the domain
				particles_remaining_v2 = particles_remaining_v2 - 1;
				total_captured_v2 = total_captured_v2 + 1;

				receptors[rec_loop].status = 0;
	        }
		}

		// update particle movement
		for(j=0;j<N;j++){

			// particle is diffusing in the domain
			if(particles[j].status==0){

				// find the next point via brownian noise
				particles[j].x_pos2 = particles[j].x_pos1+k_constant*randn(0,1);
				particles[j].y_pos2 = particles[j].y_pos1+k_constant*randn(0,1);
				particles[j].z_pos2 = particles[j].z_pos1+k_constant*randn(0,1);

				// particle hit the cylinder's side boundary and is reflected back into the domain
				if(pow(particles[j].x_pos2, 2)+pow(particles[j].y_pos2, 2)>=pow(R,2)){

					x1 = particles[j].x_pos1; y1 = particles[j].y_pos1;
					x2 = particles[j].x_pos2; y2 = particles[j].y_pos2;

					// find the intersection point on the circle
					m = (y2-y1)/(x2-x1);
					a = pow(m,2)+1;
					b = 2*m*y1-2*x1*pow(m,2);
					c = pow(m,2)*pow(x1,2)-2*m*x1*y1+pow(y1,2)-pow(R,2);

					// option one
					x_inter1 = (-b + sqrt(pow(b,2)-4*a*c))/(2*a);
					y_inter1 = m*(x_inter1-x1)+y1;

					// option two
					x_inter2 = (-b - sqrt(pow(b,2)-4*a*c))/(2*a);
					y_inter2 = m*(x_inter2-x1)+y1;

					// take the true intersection to be the closer of the two options
					d1 = pow((x2-x_inter1),2)+pow((y2-y_inter1),2);
					d2 = pow((x2-x_inter2),2)+pow((y2-y_inter2),2);
					if(d1 < d2){
						x_inter_final = x_inter1;
						y_inter_final = y_inter1;
					}else{
						x_inter_final = x_inter2;
						y_inter_final = y_inter2;
					}

					// find the m and c of y=mx+C of the tangent line
					// where the intersection occurs
					tangent_slope = -x_inter_final/y_inter_final;
					C = pow(x_inter_final,2)/y_inter_final + y_inter_final;

					// reflect the point across this tangent line
					d = (x2 + (y2 - C)*tangent_slope)/(1+pow(tangent_slope,2));
					x_reflected = 2*d-x2;
					y_reflected = 2*d*tangent_slope-y2+2*C;

					particles[j].x_pos2 = x_reflected;
					particles[j].y_pos2 = y_reflected;
				}

				/*
				** particle hit the lower boundary in the cleft
				** two options: 1) hit receptor, 2) hit reflecting region
				*/
				if(particles[j].z_pos2 <=0){

					// test to see if it hit a receptor
					for(rec_loop=0;rec_loop<M;rec_loop++)
					{
						// particle did hit a receptor and it is avaliable
				        if((pow((particles[j].x_pos2-receptors[rec_loop].center_x), 2)+
							pow((particles[j].y_pos2-receptors[rec_loop].center_y),2))<=pow(rec_rad, 2)
							&& current_time>receptors[rec_loop].switch_time)
						{

							receptors[rec_loop].switch_time = current_time + -log(((double) rand () / RAND_MAX))*tau_r;
							receptors[rec_loop].status = 1;

							// particle has been captured and is removed from the domain
							particles_remaining_v1 = particles_remaining_v1 - 1;
							total_captured_v1 = total_captured_v1 + 1;
							particles[j].status = 2;

							// exits out of the receptor loop early, since particle is now captured
							rec_loop=M+1;
				        }
					}
					// Receptor was not hit if rec_loop==M so reflect point
					if(rec_loop==M)
					{
						particles[j].z_pos2 = -particles[j].z_pos2;
					}
				}
				// particle hit upper boundary and is reflected back in
				else if(particles[j].z_pos2 >= L_z){
					particles[j].z_pos2 = 2*L_z-particles[j].z_pos2;
				}

				// update the particle's position
				particles[j].x_pos1 = particles[j].x_pos2;
				particles[j].y_pos1 = particles[j].y_pos2;
				particles[j].z_pos1 = particles[j].z_pos2;
			}

			// if all particles are removed, end trial
			if(particles_remaining_v1 == 0){
				j = N+1;
			}

		}

		// Print P(t), C(t), and R(t) every 0.0001 time units
		if(time_count == print_num)
		{
			num_recs_available=0;
			for(rec_loop=0;rec_loop<M;rec_loop++)
			{
		        if(current_time>=receptors[rec_loop].switch_time)
				{
					num_recs_available = num_recs_available+1;
		        }
			}
			
			output_matrix_v1[total_prints] = current_time;
			output_matrix_v1[total_prints+1*max_time_points] = particles_remaining_v1;
			output_matrix_v1[total_prints+2*max_time_points] = total_captured_v1;
			output_matrix_v1[total_prints+3*max_time_points] = num_recs_available;
				
			output_matrix_v2[total_prints] = current_time;			
			output_matrix_v2[total_prints+1*max_time_points] = particles_remaining_v2;
			output_matrix_v2[total_prints+2*max_time_points] = total_captured_v2;
			output_matrix_v2[total_prints+3*max_time_points] = num_recs_available;
			
			time_count =1;
			
			total_prints = total_prints + 1;
		}else
		{
			time_count = time_count + 1;
		}

		current_time +=  dt;
	}
		
	mxFree(receptors);
	mxFree(particles);
}