function [geom_out_data, out_data, dyn_data] = dmst_par_loop(sim_settings, sim_input, dmst_input)
% Single turbine plane solver. This code can be easily parallelized in a parfor loop

%% Variables declaration %%
n_ring = sim_input.n_ring;
a1 = NaN(n_ring/2,1);
a2 = a1;
delta_theta = dmst_input.delta_theta;
U_inf = dmst_input.U_inf_k;

R3 = dmst_input.R3_k';
T3 = dmst_input.T3_k';

out_data = NaN(n_ring,7);
geom_out_data = NaN(n_ring,16);
dyn_data = NaN(n_ring,8);

tsr = dmst_input.omega*unique(R3)/U_inf;

if sim_settings.st_curvature ~= 0
    sol_space = [0.5 1.2];
else
    sol_space = [0.5 1];
end
%% Evaluate first half of turbine %%
for i=1:n_ring/2
    % Define solution interval %
    sol_space_corr = [0.5*sqrt(1+(curv_corr(tsr,T3(i)))^2) sqrt(1+(curv_corr(tsr,T3(i)))^2)];
    
    % Define function to be iteratively solved %
    int_f1_eval = @(x) dmst_calc(sim_settings, sim_input, dmst_input, i, x);
    
    % Solve single upstream streamtube -> find a1(i) %
    try
        a1(i) = fzero(int_f1_eval, sol_space_corr);
    catch
        try % Increase interval range %
            a1(i) = fzero(int_f1_eval, [0 2]);
        catch % Determine interval range automatically, slower but more reliable %
            a1(i) = fzero(int_f1_eval, 1);
        end
        
        % Correct outputs if needed %
        if a1(i) > sol_space_corr(2)
            a1(i) = sol_space_corr(2);
        end
        
        if a1(i) < sol_space_corr(1)
            a1(i) = sol_space_corr(1) + 0.01;
        end
    end
    
    % Evaluate everything for a1(i) found %
    [~, geom_out_data(i,:), out_data(i,:), dyn_data(i,:)] = int_f1_eval(a1(i));
end

% If streamtube curvature correction is on, a1 needs to be recalculated %
if sim_settings.st_curvature ~= 0
    a1 = geom_out_data(1:n_ring/2,15)./U_inf;
    out_data(1:n_ring/2,1) = a1; 
end

%% Streamtube expansion %%
a1_eq = NaN(n_ring/2,1);

if sim_settings.st_expansion == 1
    % Inspired by (), basically useless %
    funz = cell(n_ring/2,1);
    funz_eq = cell(n_ring/2,1);
    funz_or = cell(n_ring/2,1);
    par_sum = NaN(n_ring/2,1);
    U_eq = NaN(n_ring/2,1);
    
    l_or = R3(1:n_ring/2) .* delta_theta .*abs(sin(T3(1:n_ring/2)));
    l_eq = R3(1:n_ring/2) .* delta_theta .*abs(sin(T3(1:n_ring/2))) ./ (2 .* a1 - 1);
    y_st = R3(1:n_ring/2) .* cos(T3(1:n_ring/2));
    
    y_eq = [y_st - l_eq/2, y_st + l_eq/2];
    y_or = [y_st - l_or/2, y_st + l_or/2];
    
    %val_x = linspace(min(y_or(:)),max(y_or(:)),1000);
    
    for i=1:n_ring/2
        funz_eq{i} = @(y) U_inf * (2 * a1(i) - 1) * ((y-y_eq(i,1)) > 0 & (y_eq(i,2)-y) > 0);
        funz{i} = @(y) 1 * ((y-y_eq(i,1)) > 0 & (y_eq(i,2)-y) > 0);
        funz_or{i} = @(y) U_inf * (2 * a1(i) - 1) * ((y-y_or(i,1)) > 0 & (y_or(i,2)-y) > 0);
    end
    %
    %     figure(2)
    %     fplot(funz_eq,[-3.5 3.5],'MeshDensity',50)
    %     figure(3)
    %     fplot(funz_or,[-3.5 3.5],':','MeshDensity',50)
    
    for i=1:n_ring/2
        %par_sum = 0;
        
        for j=1:n_ring/2
            %par_sum = par_sum + integral(funz{j},y_or(i,1),y_or(i,2));
            par_sum(j) = integral(funz{j},y_or(i,1),y_or(i,2)) / (y_or(i,2) - y_or(i,1));
        end
        
        %U_eq(i) = par_sum / (y_or(i,2) - y_or(i,1));
        U_eq(i) = sum(par_sum.*(U_inf.*(2.*a1-1)))./sum(par_sum);
        a1_eq(i) = 0.5 *(U_eq(i)/U_inf + 1);
    end
    
    % figure
    % hold on
    % grid on
    % plot(a1)
    % plot(a1_eq)
    
    % L_eq = (a1 .* R3 .* delta_theta .*abs(sin(T3-beta_k))) ./ (2 .* a1 - 1);
    % L_or = R3 .* delta_theta .*abs(sin(T3-beta_k));
    % y_st = [R3 .* cos(T3-beta_k) + L_eq/2, R3 .* cos(T3-beta_k) - L_eq/2];
    % y_or = [R3 .* cos(T3-beta_k) + L_or/2, R3 .* cos(T3-beta_k) - L_or/2];
    % v_dw = U_inf .*(2.*a1-1);
    %
    % v_dw_corr = zeros(1,n_ring);
    % for i=1:n_ring/2
    %     y_max = y_st(:,2) <= y_or(i,1);
    %     y_min = y_st(:,1) >= y_or(i,2);
    %
    %     imp_st = y_max .* y_min;
    %
    %     delta_L = imp_st.*min([y_st(:,1)-y_or(i,2) , y_or(i,1)-y_st(:,2)],[],2);
    %     delta_L(i) = L_or(i);
    %
    %     v_dw_corr(i) = sum(imp_st .* v_dw .* delta_L) ./ sum(imp_st .* delta_L);
    % end
    %
    % v_dw_corr = fliplr(v_dw_corr);
    % dmst_input.U3_k_dw = v_dw_corr .* cos(beta_k);
    % dmst_input.V3_k_dw = v_dw_corr .* sin(beta_k);
    % % modifica trascurabile... %
    
    a_upstr = flipud(a1_eq);
    
    figure(4)
    plot(y_st,U_inf.*(2.*a1-1),y_st,U_inf.*(2.*a1_eq-1),y_st,U_inf.*(2.*ones(size(a1)).*mean(a1_eq)-1));
    
elseif sim_settings.st_expansion == 2
    % Deluca, with tuning %
    l_or = R3(1:n_ring/2) .* delta_theta .*abs(sin(T3(1:n_ring/2)));
    l_eq = R3(1:n_ring/2) .* delta_theta .*abs(sin(T3(1:n_ring/2))) ./ (2 .* a1 - 1);
    y_st = R3(1:n_ring/2) .* cos(T3(1:n_ring/2));
    
    if tsr < 2.1
        f = 0.12;
    elseif tsr > 3.6
        f = 0.3;
    else
        f = -2.30308 + 1.75692*tsr - 0.287179*tsr^2;
    end
    
    %f = 1;
    
    y_eq = ((1-f)+f*(sum(l_eq)/sum(l_or))).*y_st;
    
    U_eq = interp1(y_eq,(2*a1-1).*U_inf,y_st);
    a1_eq = 0.5 * (U_eq/U_inf + 1);
    
    a_upstr = flipud(a1_eq);
    
elseif sim_settings.st_expansion == 3
    % Deluca, no tuning %
    l_or = R3(1:n_ring/2) .* delta_theta .*abs(sin(T3(1:n_ring/2)));
    l_eq = R3(1:n_ring/2) .* delta_theta .*abs(sin(T3(1:n_ring/2))) ./ (2 .* a1 - 1);
    y_st = R3(1:n_ring/2) .* cos(T3(1:n_ring/2));
    
    y_eq = (sum(l_eq)/sum(l_or)).*y_st;
    
    U_eq = interp1(y_eq,(2*a1-1).*U_inf,y_st);
    a1_eq = 0.5 * (U_eq/U_inf + 1);
    
    a_upstr = flipud(a1_eq);
    
else
    % No correction %
    a_upstr = flipud(a1);
    
end
%% Evaluate second half of turbine %%
count = 0;

for i=n_ring/2+1:n_ring
    count = count + 1;
    
    int_f2_eval = @(x) dmst_calc(sim_settings, sim_input, dmst_input, i, x, a_upstr(count));
    
    try
        a2(count) = fzero(int_f2_eval, sol_space);
    catch
        try
            a2(count) = fzero(int_f2_eval, [-5 2]);
        catch
                a2(count) = fzero(int_f2_eval, 1);
        end
        
        if a2(count) > 1
            a2(count) = 1;
        end
        
        if a2(count) < 0.5
            a2(count) = 0.5;
        end
    end
    
    [~, geom_out_data(i,:), out_data(i,:), dyn_data(i,:)] = int_f2_eval(a2(count));
end

%% Add a1_eq column to out_data %%
out_data = [out_data [a1_eq; nan(n_ring/2,1)]];

end