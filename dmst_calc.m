function [f_eff, geom_data, out_data, dyn_data] = dmst_calc(sim_settings, sim_input, dmst_input, pos_theta, a, a_upstream)
% Streamtube solver

%% Define variables %%
out_data = zeros(7,1);
geom_data = zeros(16,1);
dyn_data = NaN(8,1);

rho = sim_input.rho;
c = sim_input.c;
mu = sim_input.mu;
n_points = sim_input.n_points;
x_from_nose = sim_input.x_from_nose;
blades = sim_input.blades;

omega = dmst_input.omega;
R3 = dmst_input.R3_k(pos_theta);
T3 = dmst_input.T3_k(pos_theta);
mu_z = dmst_input.mu_z_k(pos_theta);
delta_z = dmst_input.delta_z_k;
%beta_k = deg2rad(dmst_input.beta_k);
U_inf = dmst_input.U_inf_k;

tsr = omega*R3/U_inf;

dyn_input = dmst_input.dyn_input_k;
n_ring = size(dyn_input,1);

delta_theta = 2*pi/n_ring;

%% Define new velocity %%
if nargin == 5 % Front streamtube %
    if sim_settings.st_curvature == 1
        % Deluca model, potential flow, no tuning %
        f = sin(2*T3);
        U3 = (U_inf*a)/sqrt(1+f^2);
        V3 = U3*f;
        
    elseif sim_settings.st_curvature == 2
        % Deluca model, with data from CFD %
        f = curv_corr(tsr,T3);
        
        U3 = (U_inf*a)/sqrt(1+f^2);
        V3 = U3*f;
        
    else
        % No correction %
        U3 = a*dmst_input.U3_k(pos_theta);
        V3 = a*dmst_input.V3_k(pos_theta);
        
    end
else % Back streamtube %
    U3 = a*(2*a_upstream-1)*dmst_input.U3_k(pos_theta);
    V3 = a*(2*a_upstream-1)*dmst_input.V3_k(pos_theta);
end

TSR_loc = omega * R3 / sqrt(U3^2+V3^2);

%% Geom data %%
cosTh = cos(T3);
sinTh = sin(T3);
W0 = U3*cosTh+sign(omega)*(V3*sinTh+omega*R3);
W1 = U3*sinTh-sign(omega)*V3*cosTh;
modW = sqrt(W0^2 + W1^2);
Re_c = rho*c*modW/mu;

sinAl = W1/modW;
cosAl = W0/modW;
alpha = rad2deg(atan2(W1,W0));

%% Flow curvature correction
if sim_settings.virtual == 1
    v_alpha = virtual(n_points, omega, c, x_from_nose, R3, alpha, U3, V3, sinTh, cosTh);
    alpha = alpha + v_alpha; % new, corrected alpha %
    sinAl = sin(deg2rad(alpha)); % new sinAl %
    cosAl = cos(deg2rad(alpha)); % new cosAl %
else
    v_alpha = NaN;
end

%% Lift and Drag forces assessment
dyn_prev_k = dyn_input;

% dynamic stall assessment %
if sim_settings.dynamic == 1
    % Rocchio model, old %
    [cl, alpha_d, f_d, cl_forc, cl_vortex, flag_lev, x_lev, cl_circ, ~] = ...
        cl_dyn_old_mex(modW, dyn_prev_k(pos_theta,1), dyn_prev_k(pos_theta,2), dyn_prev_k(pos_theta,3), ...
        dyn_prev_k(pos_theta,4), dyn_prev_k(pos_theta,5), dyn_prev_k(pos_theta,6), dyn_prev_k(pos_theta,7), dyn_prev_k(pos_theta,8), ...
        alpha, Re_c, n_ring, omega, c);
    
    cd = static_data_mex(alpha_d,Re_c,1); % static Cd evaluated with dynamic alpha %
    
elseif sim_settings.dynamic == 2
    % Rocchio model, new %
    alpha_in = dyn_prev_k(pos_theta,1);
    alpha_in_0 = dyn_prev_k(pos_theta,2);
    alpha_max = min([112.5*tsr^(-2.035)+3.718 50]);
    alpha_min = max([-6924*tsr^(-7.303)-4.597 -50]);
    
    cl = cl_dynamic_mex(c, omega, U_inf, modW, Re_c, alpha_in, alpha_in_0, alpha_max, alpha_min, cosTh);
    cd = static_data_mex(alpha_in,Re_c,1);
    
else
    % Static data %
    cl = static_data_mex(alpha,Re_c,0); % Cl %
    cd = static_data_mex(alpha,Re_c,1); % Cd %
    
    if sim_settings.dynamic == 1
        alpha_d = dyn_prev_k(pos_theta,1);
        f_d = dyn_prev_k(pos_theta,2);
        cl_forc = dyn_prev_k(pos_theta,3);
        cl_vortex = dyn_prev_k(pos_theta,4);
        flag_lev = dyn_prev_k(pos_theta,5);
        x_lev = dyn_prev_k(pos_theta,6);
        cl_circ = dyn_prev_k(pos_theta,7);
    end
    
end
%% Tip losses assessment
if sim_settings.tip_losses == 1
    f3 = tip_losses(blades,mu_z,sinAl); % f3 tip loss factor %
elseif sim_settings.tip_losses == 2
    f3 = tip_new_mex(mu_z,tsr,sim_input.AR);
else
    f3 = 1;
end

cl_tl = cl*f3;
cd_tl = cd*f3;

%% Streamtube Thrust coefficient evaluation %%
% Evaluated without tip losses %
if sim_settings.st_curvature ~= 0
    if nargin == 5
        a = U3/U_inf;
    end
end

ct_teo = (((5*a - 5)^2 + 729/1000)^(1/2) - 5*a + 5)^(1/3) - 9/(10*(((5*a - 5)^2 + 729/1000)^(1/2) - 5*a + 5)^(1/3)); % Spera %
%ct_teo = 4*a*(1-a); % AD theory %

ct = cl*sinAl - cd*cosAl;
cn = cl*cosAl + cd*sinAl;

cx = -(ct*cosTh - cn*sinTh);
%cy = -sign(omega)*(ct*sinTh + cn*cosTh);

%f_beta = (1/2*rho*c*modW^2)*(cx*cos(beta_k)+cy*sin(beta_k))*(delta_theta*blades/(2*pi));
f_beta = (1/2*rho*c*modW^2)*cx*(delta_theta*blades/(2*pi));

if nargin == 5 % Front streamtube %
    f_teo = ct_teo * (1/2*rho*U_inf^2) * (R3 * delta_theta * abs(sin(T3)));
    out_data(2) = NaN;
    
else % Back Streamtube %
    f_teo = ct_teo * (1/2*rho*(U_inf*(2*a_upstream-1))^2) * (R3 * delta_theta * abs(sin(T3)));
    out_data(2) = a_upstream;
end

f_eff = f_teo - f_beta;

%% Define outputs
out_data(1) = a;
%out_data(2) = a_upstream;
out_data(3) = cl_tl; % with tip losses %
out_data(4) = cd_tl; % with tip losses %
out_data(5) = f3;
out_data(6) = cl; % no tip losses %
out_data(7) = cd; % no tip losses %

geom_data(1) = R3;
geom_data(2) = cosTh;
geom_data(3) = sinTh;
geom_data(4) = W0;
geom_data(5) = W1;
geom_data(6) = modW;
geom_data(7) = Re_c;
geom_data(8) = sinAl;
geom_data(9) = cosAl;
geom_data(10) = alpha;
geom_data(11) = v_alpha;
geom_data(12) = delta_z;
geom_data(13) = T3;
geom_data(14) = TSR_loc;
geom_data(15) = U3;
geom_data(16) = V3;

if sim_settings.dynamic == 1
    dyn_data(1) = alpha_d;
    dyn_data(2) = f_d;
    dyn_data(3) = cl_forc;
    dyn_data(4) = cl_vortex;
    dyn_data(5) = flag_lev;
    dyn_data(6) = x_lev;
    dyn_data(7) = cl_circ;
    dyn_data(8) = alpha;
elseif sim_settings.dynamic == 2
    dyn_data(1) = alpha;
    dyn_data(2) = alpha_in;
    dyn_data(3) = alpha_in_0;
end

end