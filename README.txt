

To run a new set of simulations:
1) Adjust the parameters in params_DiRT_cylinder.m 
   Note: PDE-ODE parameters will updated accordingly

2) Run run_both_simulations.m and answer the prompts if any appear
   10,000 particle simulation for 1 second simulation time will take 
   ~10 minutes with 8 cores running in parallel (number of trials is 
   currently set to 8)

3) For a faster parameter search, run run_PDE_ODE_simulations.m to just 
   run the full PDE-ODE model and the MM model