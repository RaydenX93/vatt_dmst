function [vel] = init_vel(sim_input, data_geom, vel_input)
% Parses flow data for DMST simulation

%% Input read %%
nz = sim_input.nz;
n_ring = sim_input.n_ring;
H = sim_input.H;
surf_dist = sim_input.surf_dist;
r = sim_input.r;
TSR = sim_input.TSR;

Z3 = data_geom.Z3;

if nz ~= 1 % 3D simulation %
    %% Get the velocity data from CFD, previous timestep %%
    dep_data = vel_input(:,1);
    vel_u = vel_input(:,2);
    vel_v = zeros(length(dep_data),1);
    
    %% Move reference system so that (0,0,0) is in the center of the turbine and fit data %%
    dep_data2 = dep_data + (H/2 + surf_dist);
    U3 = ones(nz,n_ring).*interp1(dep_data2,vel_u,Z3(:,1),'linear','extrap');
    V3 = zeros(nz,n_ring);
    
    %% Evaluate U_inf and beta for each plane
    [U_inf_zeta, beta_zeta] = beta_calc(U3(:,1),V3(:,1));
    
    %% Evaluate undisturbed flow vel and omega %%
    U_inf = trapz(Z3(:,1),U_inf_zeta)/(max(Z3(:,1))-min(Z3(:,1)));
    [U_inf_fis, beta_fis] = infty_vel(Z3(:,1),U3(:,1),V3(:,1));
    % La potenza complessiva che investe la turbina deve
    % essere la somma della potenza da tutte le direzioni su ogni piano.
    % Perciò U_inf è la media pesata dei moduli della velocità su ogni piano
    
    %% Evaluate omega and check that it is at least 0.1 rad/s %%
    om_func = @(x) mean(x .* r ./ U_inf_zeta) - TSR;
    omega = fzero(om_func,1);
    
    TSR_plane = omega .* r ./ U_inf_zeta;
    
    vel = struct('U3',U3,'V3',V3,'U_inf_fis',U_inf_fis,'beta_fis',beta_fis,'omega',omega,'vel_u',vel_u,'vel_v',vel_v,'dep_data',dep_data,'U_inf_zeta',U_inf_zeta,'beta_zeta',beta_zeta,'TSR_plane',TSR_plane,'U_inf',U_inf);
    
else % 2D simulation %
    vel_u = vel_input;
    vel_v = 0;
    
    U3 = ones(nz,n_ring).*vel_u;
    V3 = zeros(nz,n_ring);
    
    [U_inf_zeta, beta_zeta] = beta_calc(U3(:,1),V3(:,1));
    
    U_inf = U_inf_zeta;
    beta = beta_zeta;
    
    omega = TSR * U_inf / r;
    
    vel = struct('U3',U3,'V3',V3,'U_inf',U_inf,'beta',beta,'omega',omega,'vel_u',vel_u,'vel_v',vel_v,'U_inf_zeta',U_inf_zeta,'beta_zeta',beta_zeta);
end

end