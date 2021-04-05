%%
% This file reproduces Figure 2 from 
% Handy G, Lawley SD, Revising Berg-Purcell for finite receptor kinetics,
% Biophysical Journal (2021), doi: https://doi.org/10.1016/j.bpj.2021.03.021.
%
% Loads the data from DiRT and PDE simulations via the sim_database
%
% Filename key: DiRT_cyl_N_M_R_L_tau_r.mat
%   N: number of molecules released
%   M: number of receptors/capture regions
%   R: radius of the cylinder
%   L: height/length of the cylinder
%   tau_r: receptor recharge rate
%
% sim_database{1} is the name of the DiRT data file
% sim_database{2} is the name of the PDE data file
%
% Here we consider data from N = 10000, M = 500, R = 0.1, L = 1, and a
% range of receptor recharge rates
%
% Written by Gregory Handy and Sean Lawley 04/05/2021
%%
clear; close all; clc;

%% Load the filenames of completed simulations
load('sim_database.mat')

%% plotting parameters
lw=3; % linewidth
ms=15; % marker size
sp=10; % controls marker spacing

datasets_to_plot = [8 7 6 5];
num_sets = length(datasets_to_plot);

plot_mesh=1e1;
plot_mesh2=5*1e2;

color_order = {'k','r','b','g'};

%% Make plot by looping over datasets
figure(1); clf;
for ii = 1:4
    
    % Load the data
    DiRT_data = load(sim_database{datasets_to_plot(ii),1});
    PDE_data = load(sim_database{datasets_to_plot(ii),2});
    
    % thin the data points a little for the plot
    PDE_time_thin=PDE_data.PDE_ODE_params.t(1:plot_mesh:end);
    PDE_thin = PDE_data.particles_remaining_PDE_ODE_v1(1:plot_mesh:end);
    DiRT_time_thin=DiRT_data.params.t(1:plot_mesh2:end);
    DiRT_thin=DiRT_data.particles_remaining_DiRT_v1(1:plot_mesh2:end);
    
    if ii < 3
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    
    hold all
    h(ii)=plot(PDE_time_thin,PDE_thin,'-','color',color_order{ii},'linewidth',lw);
    plot(DiRT_time_thin,DiRT_thin,'s','MarkerEdgeColor',color_order{ii},'MarkerSize',ms);%,'MarkerIndices',1:sp:length(time_thin))    
end

subplot(2,1,1)
set(gca,'fontsize',16)
xlabel('t (time) [ms]');
ylabel('Molecules remaining');
legend([h(1), h(2)],{'k_c=10^1 s^{-1}','k_c=10^2 s^{-1}'});
subplot(2,1,2)
set(gca,'fontsize',16)
xlabel('t (time) [ms]');
ylabel('Molecules remaining');
legend([h(3), h(4)],{'k_c=10^3 s^{-1}','k_c=10^4 s^{-1}'});




