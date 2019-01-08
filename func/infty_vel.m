function [U_inf, beta] = infty_vel(z, u, v)
% Evaluate averaged velocity modulus and direction %

% find average weighted u
u_mm = trapz(z(:,1),u)./(max(z)-min(z));

% find average weighted v
v_mm = trapz(z(:,1),v)./(max(z)-min(z));

% find undisturbed velocity
U_inf = sqrt((u_mm).^2+(v_mm).^2);

% find undisturbed velocity direction (beta = 0 along x direction,
% positive conter-clockwise

beta = rad2deg(wrapTo2Pi(atan2(v_mm,u_mm)));
end