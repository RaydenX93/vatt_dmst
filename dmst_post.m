function post_data = dmst_post(sim_input, out_geom_data, geom_data, out_data, data_vel)
% Prepares data for plots 

% Definitions
nz = sim_input.nz;
n_ring = sim_input.n_ring;
n_blades = sim_input.blades;

R3 = out_geom_data(:,:,1);
modW = out_geom_data(:,:,6);
sinAl = out_geom_data(:,:,8);
cosAl = out_geom_data(:,:,9);
sinTh = out_geom_data(:,:,3);
cosTh = out_geom_data(:,:,2);
delta_z = out_geom_data(:,:,12);
v_alpha = out_geom_data(:,:,11);
alpha = out_geom_data(:,:,10);

cl = out_data(:,:,3);
cd = out_data(:,:,4);
a = out_data(:,:,1);

rho = sim_input.rho;
A_ref = sim_input.A_ref;
c = sim_input.c;

omega = data_vel.omega;
U_inf = data_vel.U_inf;
U_inf_zeta = data_vel.U_inf_zeta;

%theta = rad2deg(out_geom_data(1,:,13));
theta = rad2deg(geom_data.T3(1,:));
delta_theta = 2*pi/n_ring;

%% Per evitare scoppi di Cp_2d per U_inf troppo prossimi a 0 %%
U_inf_zeta_mod = U_inf_zeta;
U_inf_zeta_picc = U_inf_zeta <= 0.1;
U_inf_zeta_mod(U_inf_zeta_picc) = 0.1;

%% Force evaluation %%
L = cl .* 0.5 .* rho .* c .* delta_z .* modW.^2;
D = cd .* 0.5 .* rho .* c .* delta_z .* modW.^2;

ct = cl.*sinAl - cd.*cosAl;
cn = cl.*cosAl + cd.*sinAl;

cx = -(ct.*cosTh - cn.*sinTh);

%% Correzione P piccole %%
%  P_picc = abs(P) < 1e-1;
%  P(P_picc) = 0;

%% Power data %%
Ft = L .* sinAl - D .* cosAl;
P = omega .* (Ft .* R3);
ref_cp = 0.5 * rho * A_ref * U_inf.^3;
%ref_cp_zeta = 0.5 .* rho .* (2 .* R3 .* delta_z) .* (ones(1,n_ring).*U_inf_zeta).^3;
ref_cp_zeta = 0.5 .* rho .* (2 .* R3 .* delta_z) .* (ones(1,n_ring).*U_inf_zeta_mod).^3;

%% One-blade Cp data per plane %%
plt_cp_theta_zeta = P ./ ref_cp_zeta;

%% One-bladed Cp data %%
if nz ~= 1
    P_vert = sum(P); % Entire turbine power along azimuthal direction %
    plt_cp_theta = P_vert./ref_cp;
else
    plt_cp_theta = P./ref_cp;
end

%% Turbine overall averaged Cp (one-blade) %%
%cp = trapz(plt_cp_theta)/n_ring;
cp_3field = repmat(plt_cp_theta,1,3);
theta_3field = [theta-360, theta, theta+360];
cp_field = interp1(theta_3field,cp_3field,[0 theta 360], 'spline');
cp = trapz(cp_field)/(n_ring);

%% Tip losses data %%
%P_plane0 = trapz(P,2)./n_ring; % Azimuthal-averaged power turbine along the blade %
P_plane0 = trapz(plt_cp_theta_zeta,2)./n_ring;
ref_tip = P_plane0(round(mean([1 nz])));

if ref_tip == 0
    ref_tip = 1;
end

if ref_tip < 0
    plt_p_mu = (-abs(P_plane0)+2*abs(ref_tip))./abs(ref_tip);
    plt_p_mu(1) = 0;
    plt_p_mu(end) = 0;
else
    plt_p_mu = P_plane0./ref_tip; % Tip losses plot %
end

plt_mu_z = geom_data.mu_z(:,1);

%% Vertically-averaged AoAs and virtualAl %%
if nz ~= 1
    plt_alpha_theta = trapz(alpha,1)./nz;
    plt_v_alpha_theta = trapz(v_alpha,1)./nz;
else
    plt_alpha_theta = alpha;
    plt_v_alpha_theta = v_alpha;
end

%% interference factor plot %%
if nz ~= 1
    plt_a_theta = trapz(a,1)./nz;
else
    plt_a_theta = a;
end

%% total turbine power and plot %%
plt_cp_tot_theta = zeros(size(plt_cp_theta));
for i = 1:n_blades
    plt_cp_tot_theta = plt_cp_tot_theta + circshift(plt_cp_theta,(i-1)*n_ring/n_blades,2);
end
cp_tot = trapz(plt_cp_tot_theta)/(n_ring);
P_tot = cp_tot * ref_cp;
% plt_cp_tot_theta = plt_cp_theta;
% cp_tot = plt_cp_theta;

%% Ulteriori conti alla Letizia %%
% Non sono sicuro che funzioni in 3D, ma serve a poco %

%diag_cp = plt_cp_theta;
if nz == 1
    diag_p = P;
    diag_a = a;
else
    diag_p = sum(P);
    diag_a = trapz(a)./nz;
end
ave_term = n_blades * delta_theta / (2*pi);
%ref_cp = 1/2*sim_input.rho*sim_input.A_ref*(data_vel.U_inf)^3;

upwind_p = sum(diag_p(1:length(diag_p)/2))*ave_term;
upwind_cp = upwind_p/ref_cp;
upwind_a = trapz(diag_a(1:length(diag_a)/2))/(length(diag_a)/2);

% downwind %
ref_cp_dw = 1/2*rho*A_ref*(U_inf*(2*upwind_a-1))^3;

downwind_p = sum(diag_p(length(diag_p)/2+1:end))*ave_term;
downwind_cp = downwind_p/ref_cp_dw;
downwind_a = trapz(diag_a(length(diag_a)/2+1:end))/(length(diag_a)/2);

cp_tot_letizia = trapz(diag_p*ave_term)/ref_cp;
cp_tot_letizia2 = sum(diag_p*ave_term)/ref_cp;

% Define output
post_data=struct( ...
    'theta',theta, ...
    'plt_cp_theta_zeta',plt_cp_theta_zeta, ...
    'plt_cp_theta',plt_cp_theta, ...
    'cp',cp, ...
    'plt_p_mu',plt_p_mu, ...
    'plt_mu_z',plt_mu_z, ...
    'plt_alpha_theta',plt_alpha_theta, ...
    'plt_v_alpha_theta',plt_v_alpha_theta, ...
    'plt_a_theta',plt_a_theta, ...
    'cp_tot',cp_tot, ...
    'plt_cp_tot',plt_cp_tot_theta, ...
    'P',P, ...
    'P_plane0',P_plane0, ...
    'P_tot',P_tot, ...
    'cx',cx, ...
    'upwind_p',upwind_p, ...
    'upwind_cp',upwind_cp, ...
    'upwind_a',upwind_a, ...
    'downwind_p',downwind_p, ...
    'downwind_cp',downwind_cp, ...
    'downwind_a',downwind_a, ...
    'cp_tot_letizia',cp_tot_letizia, ...
	'cp_tot_letizia2',cp_tot_letizia2, ...
    'U_inf_zeta_mod',U_inf_zeta_mod, ...
    'ref_cp', ref_cp, ...
    'ref_cp_zeta', ref_cp_zeta);
end