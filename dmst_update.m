function [geom_out_data, out_data, dyn_data] = dmst_update(sim_settings, sim_input, data_geom, data_vel)
% Solves one iteration of DMST on given turbine

%% Definitions %%
dynamic = sim_settings.dynamic;

R3 = data_geom.R3;
mu_z = data_geom.mu_z;
delta_z_k = data_geom.delta_z;
T3 = data_geom.T3;
delta_theta = data_geom.delta_theta;

omega = data_vel.omega;
U_inf_zeta = data_vel.U_inf_zeta;
U3 = data_vel.U3;
V3 = data_vel.V3;

nz = size(R3,1);
n_ring = size(R3,2);

geom_out_data = NaN(nz,n_ring,16);
out_data = NaN(nz,n_ring,8);
dyn_data = NaN(nz,n_ring,8);

persistent dyn_input

%% Define dyn_input %%
if isempty(dyn_input) == 1 % first iteration %
    if dynamic ~= 0
        dyn_input = zeros(nz,n_ring,8);
    else
        dyn_input = NaN(nz,n_ring,8);
    end
end

%% Run parallel solution for each turbine planes %%
parfor k=1:nz
    % Turbine "rotation" %
    if dynamic == 1
        dyn_input_k = squeeze(circshift(dyn_input(k,:,:),1,2));
    elseif dynamic == 2
        dyn_input_k = squeeze(dyn_input(k,:,:));
    else
        dyn_input_k = NaN;
    end
    
    % Define plane data %
    dmst_input = struct(...
        'R3_k',R3(k,:), ...
        'T3_k', T3(k,:), ...
        'mu_z_k', mu_z(k,:), ...
        'delta_z_k', delta_z_k, ...
        'U_inf_k', U_inf_zeta(k), ...
        'omega', omega, ...
        'U3_k', U3(k,:), ...
        'V3_k', V3(k,:), ...
        'dyn_input_k', dyn_input_k, ...
        'delta_theta', delta_theta);
    
    % Run plane solution %
    [geom_out_data(k,:,:), out_data(k,:,:), dyn_data(k,:,:)] = dmst_par_loop(sim_settings, sim_input, dmst_input);
end

% Prepare for next iteration %
if dynamic == 1
    dyn_input = dyn_data;
elseif dynamic == 2
    dyn_input = circshift(dyn_input,1,3);
    dyn_input(:,:,1) = dyn_data(:,:,1);
    dyn_input = mean(dyn_input(:,:,1:3),3).*ones(nz,n_ring,8); % to increase calculation stability %
    dyn_data = dyn_input;
end

end