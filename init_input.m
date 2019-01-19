function [sim_settings, sim_input] = init_input
% Define DMST simulation settings and turbine geometry

%% Script options %%
sim_settings = struct;
sim_settings.plot_disp = 0; % output plots display (0: text only, 1: plots) %
sim_settings.dynamic = 2; % dynamic stall model (0: static, 1: rocchio old, 2: chicchiero) %
sim_settings.tip_losses = 1; % tip losses model (0: inactive, 1: Prandtl, 2: CFD 3D) %
sim_settings.virtual = 0; % flow curvature model (0: inactive, 1: deluca) %
sim_settings.st_curvature = 2; % streamlines curvature model (0: inactive, 1: potential flow, 2: deluca, tuning) %
sim_settings.st_expansion = 2; % stramtubes expansion model (0: inactive, 1: test, 2: deluca, tuning) %
sim_settings.om_calc = 0; % omega calculation algorythm (0: normal average of TSR planes, 2: weighted average over plane power, 3: use Cp-TSR curve %
%% Turbine grid parameters %%
nz = 51; % number of cells along vertical direction (better if it is odd number) %
n_ring = 40; % number of cells along azimuthal direction (must be even number)%
                % deve essere divisibile per 2 E per il numero di pale %
n_points = 10; % number of airfoil points for flow curvature correction routine %

%% Turbine geometry and other settings (SI units)
r = 3.16228; % turbine radius %
c = 0.42162; % blade chord %
blades = 4; % n. of blades %

AR = 1.35; % aspect ratio %
TSR = 2.625; % tip-speed Ratio %
x_from_nose = c/4; % distance between airfoil leading edge and pivoting point of the profile %

surf_dist = 2; % distance between turbine top and sea surface %
bott_dist = 2; % minimum distance between turbine bottom and sea bottom %
rho = 998.2; % fluid density [kg/m3] %
mu = 0.001003; % dynamic viscosity [kg/(m s)] %

% stop criterium %
stop_it = 5; % number of iterations to consider %
stop_c = 5e-4; % stop the sim when the relative cp difference is lower than given value %

%% Evaluate quantities and produce outputs %%
% turbine height %
if nz == 1 % 2D case %
    H = 1;
else
    H = AR*(2*r);
end

A_ref = H*(2*r); % turbine frontal area %

sim_input = struct('r',r,'c',c,'blades',blades,'x_from_nose',x_from_nose,'AR',AR,'TSR',TSR,'surf_dist',surf_dist,'rho',rho,'mu',mu,'H',H,'A_ref',A_ref,'nz',nz,'n_ring',n_ring,'n_points',n_points,'stop_c',stop_c,'bott_dist',bott_dist,'stop_it',stop_it);
end