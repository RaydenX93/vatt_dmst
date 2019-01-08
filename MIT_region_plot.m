function out_fig = MIT_region_plot(lat_coord,long_coord,data,turb_mask,name)
% Creates fancy site assesment figures with m_map and cmocean %
addpath('plotting\m_map', 'plotting\cmocean');

%% Set data values on boundaries and other non-valid cells %%
% turb_mask values %
% 0: sim can be carried out %
% 1: non-physical value of velocity %
% 2: boundary %
% 3: land %
% 4: turbine intersects bottom %

data_work = data;

% Color the non-valid cells %
data_work(turb_mask == 1 | turb_mask == 2 | turb_mask == 3) = min(data(:));

% Interpolate over NaNs of valid cells %
if any(isnan(data_work(turb_mask == 0)),'all')
    good_ids = (isnan(data_work) == 0);
    long_good = long_coord(good_ids);
    lat_good = lat_coord(good_ids);
    data_good = data_work(good_ids);
    
    data_work = griddata(long_good,lat_good,data_good,long_coord,lat_coord,'natural');
    
    % On good cells where no sim has been carried out, apply min value (to
    % avoid white areas 
    data_work(isnan(data_work)) = min(data(:));
end

%% Create higher resolution grid %%
n = 1000; % Resolution %
[long_fitt, lat_fitt] = meshgrid(linspace(min(long_coord(:)),max(long_coord(:)),n),linspace(min(lat_coord(:)),max(lat_coord(:)),n));

%% Interpolate on new grid %%
data_new = griddata(long_coord,lat_coord,data_work,long_fitt,lat_fitt,'natural');

%% Initialize output figure %%
out_fig = figure('Name', name);
hold on
m_proj('miller', 'lat', [min(lat_coord(:)) max(lat_coord (:))], 'lon', [min(long_coord(:)) max(long_coord(:))]);
m_grid('box', 'fancy', 'FontSize', 15);

%% Plot data (background) %%
% Mappa dei Cp, che sarà sullo sfondo della figura %
plots = m_pcolor(long_fitt, lat_fitt, data_new);
plots.EdgeColor = 'none';
plots.FaceColor = 'interp';
colormap(cmocean('balance'));

%% Define patches where sim has not been carried out %%
turb_mask2 = turb_mask;
% The patches will be made up of areas (2), (3), (4) %
turb_mask2(turb_mask == 1) = 0; % non-valid cells "become" valid %
turb_mask2(turb_mask == 2 | turb_mask == 3 | turb_mask == 4) = 1;

B = bwboundaries(turb_mask2'); % Requires Image Processing Toolbox %

%% Plot patches %%
% Invalid cells patches %
long_vet = long_coord(:,1);
lat_vet = lat_coord(1,:);
for k = 1:length(B)
    boundary = B{k};
    m_patch(long_vet(boundary(:,2)), lat_vet(boundary(:,1))', 'black','EdgeAlpha',0);
end

% Coastline patches %
m_gshhs_f('patch', [.7 .7 .7]);

end