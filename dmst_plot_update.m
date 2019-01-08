function dmst_plot_update(fig_out, data_post)
% Creates DMST output plots
figure(fig_out)

iter = data_post.iter;
cp_iter = data_post.cp_iter;
time_iter = data_post.time_iter;

% Cp_medio-theta plot con tempo iterazioni %
sb1 = subplot(2,3,1);
yyaxis left
plt1_1 = plot(1:iter,cp_iter);
ylabel('Cp');
yyaxis right
pl11_2 = plot(1:iter,time_iter);
ylabel('Iteration time (s)');
title('One-blade turbine Cp');
xlabel('Iteration');
grid on
box on

% Interference factor plot %
sb2 = subplot(2,3,2);
plt2 = plot(data_post.theta,data_post.plt_a_theta);
ylabel('a');
title('Interference factor');
xlabel('\theta [°]');
xlim([0 360]);
xticks([0 90 180 270 360]);
grid on
box on

% One-blade turbine Cp %
sb3 = subplot(2,3,3);
plt3 = plot(data_post.theta,data_post.plt_cp_theta);
title('One-blade turbine Cp');
xlabel('\theta [°]');
ylabel('Cp');
xlim([0 360]);
xticks([0 90 180 270 360]);
grid on
box on

% Alpha-theta %
sb4 = subplot(2,3,4);
cla(sb4)
hold on
plt4_1 = plot(data_post.theta,data_post.plt_alpha_theta);
if any(isnan(data_post.plt_v_alpha_theta)) == 0
    plt4_2 = plot(data_post.theta,data_post.plt_v_alpha_theta);
end
title('Angle of attack');
xlabel('\theta [°]');
ylabel('Degrees [°]');
xlim([0 360]);
xticks([0 90 180 270 360]);
if any(isnan(data_post.plt_v_alpha_theta)) == 0
    legend('Geom. AoA','Virtual AoA');
end
grid on
box on

% Tip losses %
if isscalar(data_post.vel_input) ~= 1
    sb5 = subplot(2,3,5);
    plt5 = plot(data_post.plt_mu_z,data_post.plt_p_mu);
    title('Tip losses');
    xlabel('\mu');
    ylabel('C_P/C_p(0)');
    xlim([-1.1 1.1]);
    %ylim([0 1.15]);
    xticks([-1 -0.5 0 0.5 1]);
    grid on
end

% Velocity profile %
if isscalar(data_post.vel_input) ~= 1
    sb6 = subplot(2,3,6);
    %plt6 = plot(data_post.vel_input(:,2),data_post.vel_input(:,1));
    plt6 = plot(data_post.vel_turb(:,2),data_post.vel_turb(:,1));
    title('Incoming flow velocity profile');
    xlabel('Velocity [m/s]');
    ylabel('Depth [m]');
    grid on
end
box on

drawnow
end