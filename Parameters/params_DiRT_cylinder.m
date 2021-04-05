%%
% Non-dimensional parameters for DiRT simulations in the cylinder
%
% Written by Gregory Handy, 04/05/2021
%%
function DiRT_params = params_DiRT_cylinder()

%% Key simulation parameters that can easily be adjusted 
% Particle details
DiRT_params.n = 10000;     % number of particles
DiRT_params.D = 1;         % diffusion coefficient

% Domain dimensions
DiRT_params.L = 1.0;     % height of the cylinder
DiRT_params.half_R = 0.05;      % radius of the cylinder is 2*R
DiRT_params.rec_radius = 0.001; % radius of the receptors

% Receptor details
DiRT_params.N = 500;       % number of receptors
DiRT_params.tau_r = 1.00;  % recharge time (i.e., recharge rate is 1/tau_r)

% Numerical details
DiRT_params.num_trials = 8;
DiRT_params.max_time = 10; % max time to run the simulation,

% Define where the particles should start
% Gaussian with center at (x,y,z) and standard deviation 0.01
DiRT_params.x_start_loc = 0;
DiRT_params.y_start_loc = 0;
DiRT_params.z_start_loc = 0.9;
DiRT_params.start_std = 0.01;

%% No need to adjust the parameters below this point
% Note: simulation will automatically end once all particles have been
% removed and receptors have returned to the open state
DiRT_params.dt = 0.000001; % time step for diffusion update

DiRT_params.dt_saved = 0.001;
DiRT_params.print_num = floor(DiRT_params.dt_saved/DiRT_params.dt); % when to save the state of the system
DiRT_params.max_time_points = DiRT_params.print_num*DiRT_params.max_time;

DiRT_params.t = DiRT_params.dt_saved:DiRT_params.dt_saved:DiRT_params.max_time;
end