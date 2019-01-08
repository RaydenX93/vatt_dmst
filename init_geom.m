function [init] = init_geom(sim_input)
% Initialize turbine geometry

%% Input read & preliminary Geometry calculations %%
nz = sim_input.nz;
n_ring = sim_input.n_ring;
H = sim_input.H;
r = sim_input.r;

delta_theta = 2*pi/n_ring;
delta_z = H/nz;

R3 = ones(nz,n_ring).*r; % Straight-bladed turbine %
T3 = ones(nz,n_ring).*linspace(delta_theta/2,2*pi-delta_theta/2,n_ring);

if nz == 1
    Z3 = 0;
    mu_z = zeros(1,n_ring);
else
    Z3 = ones(nz,n_ring).*linspace(-H/2,H/2,nz)'; % Uniformly spaced turbine planes %
    mu_z = Z3./(H/2); % non-dimensional vertical position %
end

X3 = (-R3).*sin(T3);
Y3 = R3.*cos(T3);

init = struct('delta_theta',delta_theta,'delta_z',delta_z,'R3',R3,'T3',T3,'Z3',Z3,'X3',X3,'Y3',Y3,'mu_z',mu_z);
end