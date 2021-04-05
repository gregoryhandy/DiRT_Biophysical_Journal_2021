%%
% Non-dimensional parameters for PDE simulations in the cylinder
%
% Takes in the DiRT parameters and makes the necessary adjustments so that
% the two simulations match in terms of parameters
%
% Written by Gregory Handy, 04/05/2021
%%
function PDE_ODE_params = params_PDE_ODE_cylinder(params)

%% PDE-ODE parameters (uniquely determined by the DiRT parameters)

PDE_ODE_params.eps_ratio = params.rec_radius/(2*params.half_R); % ratio of receptor to cylinder radius
PDE_ODE_params.kappa = PDE_ODE_params.eps_ratio*params.N/pi;   
PDE_ODE_params.e_0 = params.N/(pi*(2*params.half_R)^2); % total receptor concentration

PDE_ODE_params.k_a = PDE_ODE_params.kappa*params.D/(PDE_ODE_params.e_0*(2*params.half_R));
PDE_ODE_params.k_b = 0;
PDE_ODE_params.k_c = 1/params.tau_r;

PDE_ODE_params.cyl_vol = (pi*(2*params.half_R)^2*params.L);

% Michaelis-Menten parameters
PDE_ODE_params.V_max = PDE_ODE_params.e_0*PDE_ODE_params.k_c/params.D;
PDE_ODE_params.K_D = (PDE_ODE_params.k_b+PDE_ODE_params.k_c)/PDE_ODE_params.k_a;

%% This parameters are the exact same as the DiRT ones
% Define where the particles should start
% Gaussian with center at z and standard deviation 0.01
PDE_ODE_params.z_start_loc=params.z_start_loc;
PDE_ODE_params.start_std = params.start_std;

% Particle details
PDE_ODE_params.n = params.n;     % number of particles
PDE_ODE_params.D = params.D;         % diffusion coefficient
PDE_ODE_params.N = params.N;       % number of receptors
PDE_ODE_params.half_R = params.half_R;      % radius of the cylinder is 2*R
PDE_ODE_params.L = params.L;     % height of the cylinder
PDE_ODE_params.tau_r = params.tau_r;


%% Numerical details for PDE-ODE system
% off centered grids
% j = 0,...,N-1 grid points, ghost points for j = -1 and N are used
PDE_ODE_params.N_z = 1000;                    % number of grid points in z
PDE_ODE_params.dz = 1/PDE_ODE_params.N_z;
PDE_ODE_params.z = (0.5*PDE_ODE_params.dz):PDE_ODE_params.dz:(-0.5*PDE_ODE_params.dz+params.L); % vector of z values, to be used for plotting

PDE_ODE_params.t = params.dt_saved:params.dt_saved:params.max_time;
end