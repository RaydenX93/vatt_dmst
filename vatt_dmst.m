function [data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(vel_input, varargin)
% VATT DMST code v1.0
% by Stefano Deluca
%
% Refer to the manual for info on how to use this code

%% Initial Cleanup %%
close all

for k = 1:length(varargin)
    if ischar(varargin{k}) 
        output_file = varargin{k};
    elseif isfloat(varargin{k})
        tsr_override = varargin{k};
    end
end

% Check output_file input %
if exist('output_file','var') == 0
    date_str = fix(clock);
    output_file = ['vatt_dmst_output_' num2str(date_str(1)) '_' num2str(date_str(2)) '_' num2str(date_str(3)) '_' num2str(date_str(4)) '_' num2str(date_str(5)) '_' num2str(date_str(6))];   
end

% Start console logging %
diary([output_file '_log.txt'])
diary on

% Add functions folder and start %
addpath('func');
disp(['VATT DMST code starting on ' char(datetime)])

%% Define simulations inputs %%
[sim_settings, sim_input] = init_input;

% Check for clearence at top and bottom %
if sim_input.nz ~= 1
    if isscalar(vel_input) == 1
        error('Scalar velocity input detected for 3D simulation')
    end
    
    if max(vel_input(:,1)) < -sim_input.surf_dist
        warning('The flow profile does not cover the rotor height. It will be extrapolated.')
    end
    
    if min(vel_input(:,1)) > -(sim_input.surf_dist + sim_input.H + sim_input.bott_dist)
        error('The turbine intersects the bottom.');
    end
else
    if isscalar(vel_input) ~= 1
        error('Velocity input must be a scalar for 2D simulation')
    end
end

% Override input TSR data %
if exist('tsr_override','var')
    sim_input.TSR_orig = sim_input.TSR;
    sim_input.TSR = tsr_override;
end

%% Define turbine geometry %%
data_geom = init_geom(sim_input);

%% Define flow field %%
data_vel = init_vel(sim_input, data_geom, vel_input);

%% Compile mex files %%
if exist(['static_data_mex.' mexext],'file') ~= 3
    cd func
    mex cl_dyn_old_mex.c
    mex cl_dynamic_mex.c
    coder -build static_data.prj
    coder -build tip_new.prj
    cd ..
end

%% Run turbine routine %%
iter = 0;
cp_iter = 0;
stopper = 1;

if sim_settings.plot_disp == 1
    fig_out = figure('Name','DMST Output');
end

while stopper > sim_input.stop_c
    % Start new iteration
    tic
    iter = iter + 1;
    
    % Solve DMST %
    [data_out_geom, data_out, data_dyn] = dmst_update(sim_settings, sim_input, data_geom, data_vel);
    
    % Evaluate outputs%
    data_post = dmst_post(sim_input, data_out_geom, data_geom, data_out, data_vel);
    
    time_iter(iter) = toc;
    cp_iter(iter) = data_post.cp;
    
    data_post.iter = iter;
    data_post.cp_iter = cp_iter;
    data_post.time_iter = time_iter;
    data_post.vel_input = vel_input;
    data_post.vel_turb = [data_geom.Z3(:,1), data_vel.U3(:,1)];
    
    % MAE stop criterion %
    if iter >= sim_input.stop_it + 1
        cp_grad = gradient(cp_iter);
        stopper = mean(abs(cp_grad(end-sim_input.stop_it:end)));
    end
    
    % Update plots %
    if sim_settings.plot_disp == 1
        dmst_plot_update(fig_out, data_post)
    end
    disp(['Iteration ' int2str(iter) ': Cp_1b = ' num2str(cp_iter(iter)) '; stopper = ' num2str(stopper,'%e') '; Time = ' num2str(time_iter(iter))]);
    
    % No while-loop if dyn stall is not active %
    if sim_settings.dynamic == 0 || iter > 30
        break;
    end
end

%% Ending tasks %%
sim_time = duration([0, 0, sum(time_iter)]);
data_post.sim_time = sim_time;

% Sort all structures in alphabetic order %
sim_settings = orderfields(sim_settings);
sim_input = orderfields(sim_input);
data_geom = orderfields(data_geom);
data_vel = orderfields(data_vel);
data_post = orderfields(data_post);

% Show dialog box on sim end %
if sim_settings.plot_disp == 1 && sim_settings.dynamic ~= 0
    msgbox(['DMST simulation completed in ' char(sim_time)])
else
    disp(['DMST simulation completed in ' char(sim_time)])
end

% Save outputs %
if sim_settings.plot_disp == 0
    fig_out = figure('Name','DMST Output');
    dmst_plot_update(fig_out, data_post)
end

savefig(fig_out,[output_file '.fig'])
save([output_file '.mat'],'data_post', 'data_geom', 'data_vel', 'data_dyn', 'sim_input', 'sim_settings');
diary off
end