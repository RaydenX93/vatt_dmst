close all
clear all
clc

%% Lancia simulazione 2D %%
%[data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(1.75);

%% Lancia simulazione 2D e specifica nome file output %%
%[data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(1.75, 'simulazione');

%% Lancia simulazione 3D e specifica nome file output %%
%depth = [0 -6 -13];
%vel_x = [0 0.25 1];
depth = linspace(0,-15);
vel_x = (depth.^2)/(1.25*100);
%depth = [0 -15];
%vel_x = [1.75 1.75];
vel_input = [depth', vel_x'];
[data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(vel_input, 'simulazione3d');

%% Lancia simulazione 3D e imponi TSR, a prescindere da init_input.m %%
%[data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(vel_input, 2, 'simulazione3d');
%[data_post, data_geom, data_vel, data_out_geom, data_out, data_dyn, sim_input, sim_settings] = vatt_dmst(vel_input, 'simulazione3d', 2);
