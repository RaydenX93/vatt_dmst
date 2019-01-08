% Process MIT .nc files prom PE CFD code%
close all
clear all
clc
addpath('func', 'plotting\cmocean');

%% NC file settings %%
% Just change the string below to the .nc file location %
pe_file = 'C:\Users\Stefano\Desktop\read_nc_data\Hz600mN010mw\pe_out.nc';

split_namefile = split(pe_file,'\');
namefile = split_namefile{end-1};

%% Gather required data from .nc file
tic
% get velocity grid
domain = getNCdata(pe_file,'vgrid3');

long_cells = size(domain,3);
lat_cells = size(domain,4);
z_cells = size(domain,2);

long_coord = squeeze(double(domain(1,1,:,:)));
lat_coord = squeeze(double(domain(2,1,:,:)));
z_coord = permute(squeeze(double(domain(3,:,:,:))),[2 3 1]);

z_maxdepth = squeeze(z_coord(:,:,end));

[lat_grid, long_grid] = meshgrid(1:lat_cells, 1:long_cells);

% get time information
sim_time = getNCdata(pe_file,'time')/3600; % hours since 2017-08-13 00:00:00 0:00, used to be seconds
sim_tstep = length(sim_time);

% get velocity data
vtot = getNCdata(pe_file,'vtot'); %hopefully, 1 is u,zonal and 2 is v,meridional

% sea/land/boundary mask %
sea_mask = getNCdata(pe_file,'landv'); % 2:sea, 1:boundary, 0:land %

plots(1) = figure;
fig(1) = pcolor(long_coord,lat_coord,sea_mask);
colormap(flipud(parula(3)))
fig(1).EdgeColor = 'none';
fig(1).FaceColor = 'texturemap';
title('Domain Mask (pcolor)')
xlabel('Longitude [°]')
ylabel('Latitude [°]')

plots(2) = figure;
fig(2) = surf(long_grid,lat_grid,sea_mask);
colormap(flipud(parula(3)))
fig(2).EdgeColor = 'none';
fig(2).Marker = 'o';
fig(2).MarkerSize = 3;
fig(2).FaceColor = 'none';
fig(2).MarkerFaceColor = 'flat';
fig(2).MarkerEdgeColor = 'black';
title('Domain Mask (surf)')
xlabel('Longitude [cell]')
ylabel('Latitude [cell]')

% Bathimetry plot %
plots(3) = figure;
fig(3) = surf(long_coord,lat_coord,z_maxdepth);
colormap(cmocean('ice'));
cb(3) = colorbar;
ylabel(cb(3),'Depth [m]');
fig(3).EdgeAlpha = 0.1;
fig(3).FaceColor = 'interp';
title('Depth Map')
xlabel('Longitude [°]')
ylabel('Latitude [°]')

toc
%% Inspect sea mask further %%
turb_mask = NaN(long_cells,lat_cells);

tic
for i = 1:long_cells
    for j = 1:lat_cells
        switch sea_mask(i,j)
            case 2 % sea %
                if any(any(abs(squeeze(vtot(1,:,i,j,:))) > 10000)) == 0
                    turb_mask(i,j) = 0; % valid cell %
                else
                    turb_mask(i,j) = 1; % non-physical value of velocity %
                end
            case 1 % boundary %
                turb_mask(i,j) = 2;
            case 0 % land %
                turb_mask(i,j) = 3;
        end
    end
end
toc

plots(4) = figure;
fig(4) = pcolor(long_coord,lat_coord,turb_mask);
colormap(parula(4))
fig(4).FaceColor = 'texturemap';
fig(4).EdgeColor = 'none';
title('Turbine Mask (pcolor)')
xlabel('Longitude [°]')
ylabel('Latitude [°]')

plots(5) = figure;
fig(5) = surf(long_grid,lat_grid,turb_mask);
colormap(parula(4))
fig(5).EdgeColor = 'none';
fig(5).Marker = 'o';
fig(5).FaceColor = 'none';
fig(5).MarkerFaceColor = 'flat';
fig(5).MarkerEdgeColor = 'black';
title('Turbine Mask (surf)')
xlabel('Longitude [cell]')
ylabel('Latitude [cell]')

%% Find average data %%
ave_data = cell(long_cells,lat_cells);
uinf_data = cell(long_cells,lat_cells);
ave_z_data = NaN(long_cells,lat_cells);

% Identify tidal period by inspecting at flow direction variation through time % 
tic
ave_periods = NaN(long_cells,lat_cells);
for i = 1:long_cells
    for j = 1:lat_cells
        if turb_mask(i,j) == 0
            [~ , angolo] = infty_vel(squeeze(z_coord(i,j,:)), squeeze(vtot(1,:,i,j,:)), squeeze(vtot(2,:,i,j,:)));
            [~,locs] = findpeaks(angolo);
            ave_periods(i,j) = mean(diff(locs));
        end
    end
end
tidal_period = round(mean(ave_periods,'all','omitnan')) + 1; % +1 is a correction %
toc

% Average flow data over last tidal period, considering flow modulus only %
tic
for i = 1:long_cells
    for j = 1:lat_cells
        if turb_mask(i,j) == 0         
            ave_data{i,j}(:,1) = squeeze(z_coord(i,j,:));
            
            uinf_data{i,j} = squeeze(sqrt(double(vtot(1,:,i,j,end-(tidal_period-1):end)).^2+double(vtot(2,:,i,j,end-(tidal_period-1):end)).^2)/100);
            ave_data{i,j}(:,2) = trapz(uinf_data{i,j},2)/tidal_period;
            
            ave_z_data(i,j) = trapz(-ave_data{i,j}(:,1),ave_data{i,j}(:,2))/(max(ave_data{i,j}(:,1))-min(ave_data{i,j}(:,1)));
            
        else
            ave_data{i,j} = NaN;
            uinf_data{i,j} = NaN;
        end
    end
end
toc

plots(6) = MIT_region_plot(lat_coord,long_coord,ave_z_data,turb_mask,'');
cb(6) = colorbar('FontSize', 15);
ylabel(cb(6),'Velocity [m/s]', 'FontSize', 15);
title('Averaged undisturbed velocity')

%% Save variables to file %%
disp('Saving data. This may take a few minutes')
tic

savefig(plots,[namefile '_processed_figures.fig'],'compact');

clearvars i j plots fig cb
save([namefile '_rawdata.mat'], 'vtot', 'domain','-v7.3');
clear vtot domain

save([namefile '_processed.mat'],'-v7.3');
toc

function data = getNCdata(file,var)
% Function to retrieve 'var' variable fron generic 'file' .nc
    ncid = netcdf.open(file);
    varid = netcdf.inqVarID(ncid,var);
    data = netcdf.getVar(ncid,varid);
    netcdf.close(ncid)
end