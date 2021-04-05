%%
% Note: This file only runs the PDE-ODE and MM simulations
%
% To run a new sim, adjust the parameters in the params_DiRT_cylinder.m
% file. The PDE-ODE parameters are automatically updated. 
%
% Plots and compares version 1 (particles removed from the domain when 
% they hit a receptor) with version 2 (particles removed after the receptor
% recharges)
%
% The sim_database will not be updated, as this requires both the DiRT
% simulation and PDE-ODE simulation to be completed
%
% Plots version 1 (particles removed from the domain when they hit a 
% receptor)
% Note: version 2 (particles removed after the receptor recharges) also 
% avaliable for plotting
%
% Written by Gregory Handy, 04/05/2021
%%
clear; close all; clc;

% Add all subfolders to path
% Determine where your m-file's folder is.
restoredefaultpath; 
folder = fileparts(which('run_PDE_ODE_simulations.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(folder)

%% Load parameters
DiRT_params = params_DiRT_cylinder();
PDE_ODE_params = params_PDE_ODE_cylinder(DiRT_params);

%%

data_dir = './PDE_ODE_Data/';
temp_name = sprintf('PDE_ODE_cyl_%d_%d_%.2f_%.2f_%.2f.mat',DiRT_params.n,DiRT_params.N,...
    2*DiRT_params.half_R,DiRT_params.L,DiRT_params.tau_r);

PDE_ODE_data=[];
x = [];
try
    PDE_ODE_data = load(strcat(data_dir,temp_name));
catch
    fprintf('No PDE-ODE data found for those parameter values\n');
    fprintf('Running simulations (this may take a moment) \n');
    x = 'y';
end
    
if ~isempty(PDE_ODE_data)
   fprintf('PDE-ODE data already found for that parameter set \n');
   prompt = 'Rerun simulations? (y/n) ';
   
   while ~(strcmp(x,'y') || strcmp(x,'n'))
       x = input(prompt,'s');
   end
end

if strcmp(x,'y')
    PDE_ODE_sim(PDE_ODE_params);
    PDE_ODE_data = load(strcat(data_dir,temp_name));
end

%% Plot the results

color_scheme = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];

figure(1); clf;
hh(1) = plot(PDE_ODE_params.t,PDE_ODE_data.particles_remaining_PDE_ODE_v1,'linewidth',1.5,'color',color_scheme(1,:));
hold on
hh(2) = plot(PDE_ODE_params.t,PDE_ODE_data.particles_remaining_MM,'linewidth',1.5,'color',color_scheme(2,:));

xlabel('Time (A.U.)')
ylabel('Number of Particles Remaining')
set(gca,'fontsize',16)
legend(hh, {'PDE-ODE','MM Model'})





