clear all
close all
clc

[x,fval,exitflag,output] = fminbnd(vatt_dmst_optim,x1,x2);

 function cp = vatt_dmst_optim(x)
    [data_post, ~, ~, ~, ~, ~, ~, ~] = vatt_dmst(vel_input, 'simulazione3d_optim', x);
    cp = data_post.cp_tot;
 end