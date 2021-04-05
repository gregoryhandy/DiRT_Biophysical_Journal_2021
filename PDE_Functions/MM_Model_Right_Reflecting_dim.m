%%
% Runs the continuum PDE model with Michaleis-Menten boundary conditions
%
% See the following paper for more details:
% Handy G, Lawley SD, Revising Berg-Purcell for finite receptor kinetics,
% Biophysical Journal (2021), doi: https://doi.org/10.1016/j.bpj.2021.03.021.
%
% Written by Gregory Handy, 04/05/2021
%%
function dudt = MM_Model_Right_Reflecting_dim(t, u0, V_max, K_D, D, N_x, h)

dudt = zeros(N_x+1,1);

%% ODES 

% account for the nonlinear term for the left endpoint
% technically we have two roots, but the + is clearly the right one
% test: take u0(1) = 0;
g_u = (-(2*K_D+V_max*h) + sqrt((2*K_D+h*V_max)^2-4*(-2*K_D*u0(1) ...
    -u0(1)^2 + h*V_max*u0(1))))/2;


% ODE for the left end point (partially absorbing)
dudt(1) = D/h^2*(g_u-2*u0(1)+u0(2));

% ODEs for the interior points
for i = 2:(N_x-1)
    dudt(i) = D/h^2*(u0(i-1) - 2*u0(i) + u0(i+1));
end

% ODE for the right end point (reflecting)
dudt(N_x) = D/h^2*(u0(N_x-1) - u0(N_x));

% this is the ODE for the number captured (auxiliary variable/counter)
dudt(N_x+1) = D*(u0(1)-g_u)/h;

end