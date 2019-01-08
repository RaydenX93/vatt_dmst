% Launch MIT simulation %
close all
clear all
clc
addpath('func', 'plotting\m_map', 'plotting\cmocean');

%% Load MIT processed data %%
nome_file = 'Hz600mN010mw_processed'; % Insert processed MIT data here, without extension %
dirname = [nome_file '_outputs'];
load([nome_file '.mat'])

if exist(dirname,'dir') ~= 7
    mkdir(dirname)
end

%% Initialize sim %%
sim_step = 5; % skip factor over grid -> 1 simulation every X cells %

[~, sim_input] = init_input;

dmst_run = tic;
cp = NaN(long_cells,lat_cells);
p = NaN(long_cells,lat_cells);

% Check for intersection with sea bottom %
intersec_ind = (z_maxdepth > -(sim_input.surf_dist + sim_input.H + sim_input.bott_dist)) & turb_mask == 0;
turb_mask(intersec_ind) = 4; % Add new value to turb_mask where bottom is intersected %

points_total = sum(sum(turb_mask == 0)); % number of possible simulations in the domain %
good_cells = [];

% Define grid %
for_x = 1:sim_step:size(long_coord,1);
for_y = 1:sim_step:size(lat_coord,2);

% Add last row/column in the grid if missing %
if sum(for_x == size(long_coord,1)) == 0
    for_x = [for_x size(long_coord,1)];
end

if sum(for_y == size(lat_coord,2)) == 0
    for_y = [for_y size(lat_coord,2)];
end

% Define cells where simulation can be run %
sub_mask = zeros(long_cells,lat_cells);

% Mask for subsampled grid
for i = for_x
    for j = for_y
        sub_mask(i,j) = 1;
    end
end
sub_mask = logical(sub_mask);

% Mask for areas at high averaged velocity %
vel_mask = ave_z_data >= (2/3)*max(ave_z_data,[],'all');

% Define mask on which to run sims %
run_mask = (sub_mask | vel_mask) & turb_mask == 0;

for i = 1:long_cells
    for j = 1:lat_cells
        if run_mask(i,j) == 1
            good_cells = [good_cells; i, j]; % indices of cells where simulation can be run %
        end
    end
end

%run_mask = logical(run_mask);
expected_sims = length(good_cells);
run_number = 0;

% Sub-sampling grid plot %
figure('Name', 'Sub-sampled grid');
hold on
m_proj('miller', 'lat', [min(lat_coord(:)) max(lat_coord (:))], 'lon', [min(long_coord(:)) max(long_coord(:))]);
m_grid('box', 'fancy', 'FontSize', 15);
plt = m_pcolor(long_coord, lat_coord, turb_mask);
colormap(parula(length(unique(turb_mask))))
plt.FaceColor = 'texturemap';
plt.EdgeColor = 'none';
m_scatter(long_coord(run_mask),lat_coord(run_mask));
cb(1) = colorbar('FontSize', 15);

disp('Press any key to continue')
pause
close all

%% Run simulations on grid %%
time_run = NaN(length(expected_sims));
tsr_opt = NaN(long_cells,lat_cells);
exitflag = NaN(long_cells,lat_cells);

for i=1:expected_sims
    time_run0 = tic;
    run_number = run_number + 1;
    lon = good_cells(i,1);
    lat = good_cells(i,2);
    output_name = [nome_file '_sub_' num2str(lon) '_' num2str(lat)];
    
    vel_input = [ave_data{lon,lat}(:,1), ave_data{lon,lat}(:,2)];
    
    if exist([dirname '\' output_name '.mat'],'file') ~= 0
        load([dirname '\' output_name '.mat'])
        
        if isfield(sim_input,'TSR_orig')
            exitflag(i,j) = 1;
        else
            exitflag(i,j) = 0;
        end
        
    else
        clc
        disp(['DMST Simulation running: ' num2str(i) '/' num2str(expected_sims) ' (' num2str(round((i-1)*100/expected_sims)) '%) ETA: ' char(duration([0, 0, (expected_sims-(i-1))*mean(time_run,'omitnan')]))])
        disp(['Longitude: ' num2str(lon) '; Latitude: ' num2str(lat)])
        
        % Optimize TSR
        disp('Starting TSR optimization')
        fun = @(x) vatt_dmst_optim(vel_input,x);
        [tsr_opt(i,j), ~, exitflag(i,j)] = fminbnd(fun,2,3); % TSR_min = 2, TSR_max = 3 %
        
        if exitflag(i,j) == 1
            disp('Optimization succesfully completed.')
            [data_post, ~, ~, ~, ~, ~, ~, ~] = vatt_dmst(vel_input, [dirname '\' output_name], tsr_opt(i,j));
        else
            disp('Optimization failed. Using default TSR...')
            [data_post, ~, ~, ~, ~, ~, ~, ~] = vatt_dmst(vel_input, [dirname '\' output_name]); % Use default TSR %
        end
        
    end
    
    cp(lon,lat) = data_post.cp_tot;
    p(lon,lat) = data_post.cp_tot * data_post.ref_cp;
    
    time_run(i) = toc(time_run0);
end

time_end = toc(dmst_run);
dmst_time = duration([0, 0, time_end]);

%% Make output plot %%
close all

out_fig(1) = figure('Name', 'Sub-sampled grid'); % Sub-sampling grid plot %
hold on
m_proj('miller', 'lat', [min(lat_coord(:)) max(lat_coord (:))], 'lon', [min(long_coord(:)) max(long_coord(:))]);
m_grid('box', 'fancy', 'FontSize', 15);
plt = m_pcolor(long_coord, lat_coord, turb_mask);
colormap(parula(length(unique(turb_mask))))
plt.FaceColor = 'texturemap';
plt.EdgeColor = 'none';
m_scatter(long_coord(run_mask),lat_coord(run_mask));
cb(1) = colorbar('FontSize', 15);
%ylabel(cb(1),'Turbine Mask', 'FontSize', 15);

out_fig(2) = MIT_region_plot(lat_coord,long_coord,cp,turb_mask,'Cp'); % Cp plot %
cb(2) = colorbar('FontSize', 15);
%ylabel(cb(2),'C_P', 'FontSize', 15);

out_fig(3) = MIT_region_plot(lat_coord,long_coord,p/1000,turb_mask,'P'); % P plot in kW %
cb(3) = colorbar('FontSize', 15);
%ylabel(cb(3),'P (kW)', 'FontSize', 15);

%% Save output %%
save([dirname '\' output_name '_sub_finished.mat'],'long_coord','lat_coord','cp','p','dmst_time','run_number','good_cells','turb_mask','run_mask', 'tsr_opt', 'exitflag');
save([dirname '\' output_name '_sub_finished_paranoid.mat'])

for kk = 1:3
    FigHandle = out_fig(kk);
    FigName   = get(FigHandle, 'Name');
    savefig(FigHandle, [dirname '\' output_name '_sub_finished_' FigName '.fig'],'compact');
end

disp(['Finished on ' char(datetime('now')) ' in ' char(dmst_time) ' (' char(dmst_time/run_number) '/sim)!'])

function cp = vatt_dmst_optim(vel_input, x)
[data_post, ~, ~, ~, ~, ~, ~, ~] = vatt_dmst(vel_input, 'simulazione3d_optim', x);
cp = - data_post.cp_tot; % negative because fminbnd looks for a minimum %
end