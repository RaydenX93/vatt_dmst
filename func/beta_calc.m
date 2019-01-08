function [vel_inf, beta] = beta_calc(u_mm,v_mm)
% Evaluates vector of velocity modules and directions %

% find undisturbed velocity
vel_inf = sqrt((u_mm).^2+(v_mm).^2);

% find undisturbed velocity direction (beta = 0 along x direction,
% positive conter-clockwise
beta = rad2deg(wrapTo2Pi(atan2(v_mm,u_mm)));
end