%%
% Sets up the PDE simulation, which uses ode15s 
%
% See Full_Model_Right_Reflecting_dim.m in the PDE_Functions folder for the
% ODEs being solved
%
% Written by Gregory Handy, 04/05/2021
%%
function PDE_ODE_sim(PDE_ODE_params)

%% Specify initial conditions (Guassian)
u0 = zeros(PDE_ODE_params.N_z+1,1); 
mu = PDE_ODE_params.z_start_loc; sigma_std = PDE_ODE_params.start_std;
for i=1:PDE_ODE_params.N_z   
   u0(i) = exp(-(PDE_ODE_params.z(i)-mu)^2/(2*sigma_std^2)) ...
       / sqrt(2*pi*sigma_std^2)*PDE_ODE_params.n/PDE_ODE_params.cyl_vol;
end
u0(PDE_ODE_params.N_z+1) = 0;
u0(PDE_ODE_params.N_z+2) = 0;

%% Run the full model and MM approximation

fprintf('Running PDE-ODE system \n')
options = odeset('AbsTol', 10^-6, 'RelTol', 10^-6, 'MaxStep', 0.01);
[~,u_PDE_ODE_ODE45] = ode15s(@Full_Model_Right_Reflecting_dim, PDE_ODE_params.t, u0, ...
    options, PDE_ODE_params.k_a, PDE_ODE_params.k_b, PDE_ODE_params.k_c, PDE_ODE_params.D, PDE_ODE_params.e_0, ...
    PDE_ODE_params.N_z, PDE_ODE_params.dz);

particles_remaining_PDE_ODE_v1 = PDE_ODE_params.n-u_PDE_ODE_ODE45(:,PDE_ODE_params.N_z+2)*PDE_ODE_params.cyl_vol;
particles_remaining_PDE_ODE_v2 = (PDE_ODE_params.n-u_PDE_ODE_ODE45(:,PDE_ODE_params.N_z+2)*PDE_ODE_params.cyl_vol...
    +u_PDE_ODE_ODE45(:,PDE_ODE_params.N_z+1)*PDE_ODE_params.N/PDE_ODE_params.e_0);

frac_open_receptors = u_PDE_ODE_ODE45(:,PDE_ODE_params.N_z+1)*PDE_ODE_params.N/PDE_ODE_params.e_0;


fprintf('Running Michaelis-Menten approximation \n')
[~,u_MM_ODE45] = ode15s(@MM_Model_Right_Reflecting_dim, PDE_ODE_params.t, u0(1:(PDE_ODE_params.N_z+1)), ...
    options, PDE_ODE_params.V_max, PDE_ODE_params.K_D, PDE_ODE_params.D, PDE_ODE_params.N_z, PDE_ODE_params.dz);
particles_remaining_MM = PDE_ODE_params.n-u_MM_ODE45(:,PDE_ODE_params.N_z+1)*PDE_ODE_params.cyl_vol;

%% Save the data

fprintf('Saving the data \n')
mkdir('PDE_ODE_Data');

temp_filename = sprintf('./PDE_ODE_Data/PDE_ODE_cyl_%d_%d_%.2f_%.2f_%.2f.mat',...
    PDE_ODE_params.n,PDE_ODE_params.N,2*PDE_ODE_params.half_R,PDE_ODE_params.L,PDE_ODE_params.tau_r);

save(temp_filename,'particles_remaining_PDE_ODE_v1','particles_remaining_PDE_ODE_v2',...
    'particles_remaining_MM','frac_open_receptors','PDE_ODE_params','u_MM_ODE45')

end

