%%
% Runs the DiRT simulations using mex
%
% Setup to run the trials in parallel
%
% Written by Gregory Handy, 04/05/2021
%%
function DiRT_sim_mex(params)

%% Compile the mex code
fprintf('Compiling mex code \n');
mex DiRT_cylinder_mex.c

%% Preallocate memory 
output_v1 = zeros(params.max_time_points,4,params.num_trials);
output_v2 = zeros(params.max_time_points,4,params.num_trials);
%% Loop over trials
% If unable to run in parallel, simply change this to a for-loop
tic;
fprintf('Running DiRT trials \n');
parfor hh = 1:params.num_trials
    
    %% Randomly distribute the receptors on the bottom of the cylinder
    [rec_centers] = place_receptors(2*params.half_R,params.N,params.rec_radius,0);
    
    %% Run the simulation
    [output_v1(:,:,hh), output_v2(:,:,hh)] = DiRT_cylinder_mex(params.max_time,...
        params.dt, params.D, 2*params.half_R, params.L, params.x_start_loc,...
        params.y_start_loc, params.z_start_loc, params.n, params.tau_r, ...
        params.N, rec_centers, params.rec_radius, ...
        params.start_std, params.print_num, params.max_time_points);
end
toc;

%% Adjust the output length to be consistent across trials
% Some trials might end before max_time is reached (i.e. all particles have 
% escaped and all receptors recharged). This repeats the end state for the
% remaining time bins
for hh = 1:params.num_trials
 
    temp_index = find(output_v1(:,1,hh)==0,1);
    if ~isempty(temp_index)
        remaining_slots = length(params.t)-temp_index+1;
        output_v1(temp_index:end,2:end,hh) = repmat(output_v1(temp_index-1,2:end,hh),remaining_slots,1);
    end
    
    % fixes a numerical error in the time recording
    output_v1(:,1,hh) = params.t;    
    
    temp_index = find(output_v2(:,1,hh)==0,1);
    if ~isempty(temp_index)
        remaining_slots = length(params.t)-temp_index+1;
        output_v2(temp_index:end,2:end,hh) = repmat(output_v2(temp_index-1,2:end,hh),remaining_slots,1);
    end
    
    % fixes a numerical error in the time recording
    output_v2(:,1,hh) = params.t; 
end
%%

particles_remaining_DiRT_v1 = mean(output_v1(:,2,:),3);
open_receptors_DiRT_v1 = mean(output_v1(:,4,:),3);
standard_error_v1 = sqrt(var(output_v1(:,2,:),0,3))/sqrt(params.num_trials);

particles_remaining_DiRT_v2 = mean(output_v2(:,2,:),3);
open_receptors_DiRT_v2 = mean(output_v2(:,4,:),3);
standard_error_v2 = sqrt(var(output_v2(:,2,:),0,3))/sqrt(params.num_trials);


%% Save the results
fprintf('Saving results \n');
mkdir('DiRT_Data');
temp_name = sprintf('./DiRT_Data/DiRT_cyl_%d_%d_%.2f_%.2f_%.2f.mat',params.n,params.N,...
    2*params.half_R,params.L,params.tau_r);
save(temp_name,'particles_remaining_DiRT_v1','open_receptors_DiRT_v1',...
    'standard_error_v1','particles_remaining_DiRT_v2','open_receptors_DiRT_v2',...
    'standard_error_v2','params');
end

