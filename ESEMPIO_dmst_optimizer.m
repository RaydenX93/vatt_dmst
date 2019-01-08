clear all
close all
clc

depth = [0 -15];
vel_x = [1.75 1.75];
vel_input = [depth', vel_x'];

fun = @(x) vatt_dmst_optim(vel_input,x);
[tsr_opt,fval,exitflag,output] = fminbnd(fun,2,3);

 function cp = vatt_dmst_optim(vel_input, x)
    [data_post, ~, ~, ~, ~, ~, ~, ~] = vatt_dmst(vel_input, 'simulazione3d_optim', x);
    cp = - data_post.cp_tot; % negative because fminbnd looks for a minimum %
 end
