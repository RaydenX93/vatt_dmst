clear all
close all
clc

depth = [0 -15];
vel_x = [1.75 1.75];
vel_input = [depth', vel_x'];

min_tsr = 2;
max_tsr = 3;

%% Identifica TSR ottimale per simulazione 3D %%
fun = @(x) vatt_dmst_optim(vel_input,x);
[tsr_opt,fval,exitflag,output] = fminbnd(fun,min_tsr,max_tsr); 

tsr_opt
fval = - fval

 function cp = vatt_dmst_optim(vel_input, x)
    [data_post, ~, ~, ~, ~, ~, ~, ~] = vatt_dmst(vel_input, 'simulazione3d_optim', x);
    cp = - data_post.cp_tot; % negative because fminbnd looks for a minimum %
 end
