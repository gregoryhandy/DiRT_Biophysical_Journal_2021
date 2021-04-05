%%
% This file runs the DiRT, full PDE-ODE, and MM simulations for the given
% parameter set found in Parameters folder
%
% Model details found in:
% Handy G, Lawley SD, Revising Berg-Purcell for finite receptor kinetics,
% Biophysical Journal (2021), doi: https://doi.org/10.1016/j.bpj.2021.03.021.
%
% To run a new sim, adjust the parameters in the params_DiRT_cylinder.m
% file. The PDE-ODE parameters are automatically updated. 
%
% The script checks to see if sims already exist for the chosen parameter
% set, 
%   If it does, it will asks if you want to rerun the sims
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
folder = fileparts(which('run_DiRT_and_PDE_simulations.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(folder)

%% Load parameters
DiRT_params = params_DiRT_cylinder();
PDE_ODE_params = params_PDE_ODE_cylinder(DiRT_params);

%% Check to see if DiRT data already exists
data_dir = './DiRT_Data/';
temp_name = sprintf('DiRT_cyl_%d_%d_%.2f_%.2f_%.2f.mat',DiRT_params.n,DiRT_params.N,...
    2*DiRT_params.half_R,DiRT_params.L,DiRT_params.tau_r);

DiRT_data=[];
x = [];
try
    DiRT_data = load(strcat(data_dir,temp_name));
catch
    fprintf('No DiRT data found for those parameter values\n');
    fprintf('Running DiRT simulations (this may take a moment) \n');
    x = 'y';
end

if ~isempty(DiRT_data)
   fprintf('DiRT data already found for that parameter set \n');
   prompt = 'Rerun simulations? (y/n) ';
   
   while ~(strcmp(x,'y') || strcmp(x,'n'))
       x = input(prompt,'s');
   end
end

if strcmp(x,'y')
    DiRT_sim_mex(DiRT_params);
    DiRT_data = load(strcat(data_dir,temp_name));
end

%% Check to see if the PDE data already exists

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

%% Update the simulation database
update_sim_database();

%% Plot the results

color_scheme = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];

figure(1); clf;
hh(1) = plot(PDE_ODE_params.t,PDE_ODE_data.particles_remaining_PDE_ODE_v1,'linewidth',1.5,'color',color_scheme(1,:));
hold on
hh(2) = plot(PDE_ODE_params.t,PDE_ODE_data.particles_remaining_MM,'linewidth',1.5,'color',color_scheme(2,:));
hh(3) = plot(DiRT_params.t,DiRT_data.particles_remaining_DiRT_v1,'linewidth',1.5,'color',color_scheme(3,:));

x_temp = PDE_ODE_params.t;
curve1 = (DiRT_data.particles_remaining_DiRT_v1 + DiRT_data.standard_error_v1);
curve2 = (DiRT_data.particles_remaining_DiRT_v1 - DiRT_data.standard_error_v1);
x2 = [x_temp, fliplr(x_temp)];
inBetween = [curve1', fliplr(curve2')];
fill(x2, inBetween, color_scheme(3,:),'FaceAlpha',0.2,'EdgeColor',color_scheme(3,:));
xlabel('Time')
ylabel('Number of Particles Remaining')
set(gca,'fontsize',16)
legend(hh, {'PDE-ODE','MM Model','DiRT'})






