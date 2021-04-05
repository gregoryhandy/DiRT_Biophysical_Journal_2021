%%
% Runs the continuum PDE model with boundary conditions that 
% depend on an ODE
%
% See the following paper for more details:
% Handy G, Lawley SD, Revising Berg-Purcell for finite receptor kinetics,
% Biophysical Journal (2021), doi: https://doi.org/10.1016/j.bpj.2021.03.021.
%
% Written by Gregory Handy, 04/05/2021
%%
% function dudt = Full_Model_Right_Reflecting_dim(t, u0, params, num_params)
function dudt = Full_Model_Right_Reflecting_dim(t, u0, k_a, k_b, k_c, D, e_0, N_x, h)

dudt = zeros(N_x+2,1);
 
%% ODES 
c = u0(N_x+1);
% f_w = (1-gamma*h/2*(1-w))/(1+gamma*h/2*(1-w));
f_c = (1-h*k_a/(2*D)*(e_0-c))/(1+h*k_a/(2*D)*(e_0-c));

% ODE for the left end point
dudt(1) = D/h^2*((f_c -2)*u0(1) + u0(2));

% ODEs for the interior points
for i = 2:(N_x-1)
    dudt(i) = D/h^2*(u0(i-1) - 2*u0(i) + u0(i+1));
end

% ODE for the right end point
dudt(N_x) = D/h^2*(u0(N_x-1) - u0(N_x));

% this is the ODE for c(t)
dudt(N_x+1) =  k_a*(e_0-c)*u0(1)/2*(f_c+1)-(k_b+k_c)*c;

% this is the ODE for the number captured (auxiliary variable/counter)
dudt(N_x+2) = D*u0(1)*(1-f_c)/h;

end