function [vel] = init_vel(sim_input, sim_settings,data_geom, vel_input)
% Parses flow data for DMST simulation

%% Input read %%
nz = sim_input.nz;
n_ring = sim_input.n_ring;
H = sim_input.H;
surf_dist = sim_input.surf_dist;
r = sim_input.r;
TSR = sim_input.TSR;
rho = sim_input.rho;

Z3 = data_geom.Z3;
delta_z = data_geom.delta_z;

if nz ~= 1 % 3D simulation %
    %% Get the velocity data from CFD, previous timestep %%
    dep_data = vel_input(:,1);
    vel_u = vel_input(:,2);
    vel_v = zeros(length(dep_data),1);
    
    %% Move reference system so that (0,0,0) is in the center of the turbine and fit data %%
    dep_data2 = dep_data + (H/2 + surf_dist);
    U3 = ones(nz,n_ring).*interp1(dep_data2,vel_u,Z3(:,1),'linear','extrap');
    V3 = zeros(nz,n_ring);
    
    %% Evaluate U_inf and beta for each plane
    [U_inf_zeta, beta_zeta] = beta_calc(U3(:,1),V3(:,1));
    
    %% Evaluate undisturbed flow vel and omega %%
    U_inf = trapz(Z3(:,1),U_inf_zeta)/(max(Z3(:,1))-min(Z3(:,1)));
    [U_inf_fis, beta_fis] = infty_vel(Z3(:,1),U3(:,1),V3(:,1));
    % La potenza complessiva che investe la turbina deve
    % essere la somma della potenza da tutte le direzioni su ogni piano.
    % Perciò U_inf è la media pesata dei moduli della velocità su ogni piano
    
    %% Evaluate omega and check that it is at least 0.1 rad/s %%
    %om_func = @(x) mean(x .* r ./ U_inf_zeta) - TSR;
    %omega = fzero(om_func,1);
    
    if sim_settings.om_calc == 1 % weighed average over plane power %
        omega = (TSR / r) * (sum(U_inf_zeta.^3) / sum(U_inf_zeta.^2));
    elseif sim_settings.om_calc == 2 % use Cp-TSR curve %
        disp('Optimizing omega value...')
        cfd2d_cp = [0.363; 0.387; 0.413; 0.421; 0.426; 0.42; 0.409; 0.395; 0.38; 0.366; 0.284];
        cfd2d_tsr = [2.1; 2.2; 2.3; 2.35; 2.5; 2.625; 2.75; 2.875; 3; 3.1; 3.6];
        p_teo =@(x) -1*turb_pot_cfd2d(x,delta_z,r,U_inf_zeta,rho,cfd2d_cp,cfd2d_tsr); % minus is to convert maximum problem into minimum problem %
        
        rpm_min = 1;
        rpm_max = 10;
        
        [omega, ~, exitflag] = fminbnd(p_teo,rpm_min*(2*pi/60),rpm_max*(2*pi/60));
        
        if exitflag ~= 1
            error('Optimal omega could not be found.')
        end
        
    else % normal average of TSR planes %
        omega = (TSR / r) * (nz / sum(1./U_inf_zeta));
    end   
    
    TSR_plane = omega .* r ./ U_inf_zeta;
    
    vel = struct('U3',U3,'V3',V3,'U_inf_fis',U_inf_fis,'beta_fis',beta_fis,'omega',omega,'vel_u',vel_u,'vel_v',vel_v,'dep_data',dep_data,'U_inf_zeta',U_inf_zeta,'beta_zeta',beta_zeta,'TSR_plane',TSR_plane,'U_inf',U_inf);
    
else % 2D simulation %
    vel_u = vel_input;
    vel_v = 0;
    
    U3 = ones(nz,n_ring).*vel_u;
    V3 = zeros(nz,n_ring);
    
    [U_inf_zeta, beta_zeta] = beta_calc(U3(:,1),V3(:,1));
    
    U_inf = U_inf_zeta;
    beta = beta_zeta;
    
    omega = TSR * U_inf / r;
    
    vel = struct('U3',U3,'V3',V3,'U_inf',U_inf,'beta',beta,'omega',omega,'vel_u',vel_u,'vel_v',vel_v,'U_inf_zeta',U_inf_zeta,'beta_zeta',beta_zeta);
end

end

function ptot = turb_pot_cfd2d(omega,delta_z,r,U_inf_zeta,rho,cfd2d_cp,cfd2d_tsr)
    % Evaluates the turbine power by assuming a cfd 2d Cp on every turbine
    % plane
    tsr = omega.*r./U_inf_zeta;
    cp_z = interp1(cfd2d_tsr,cfd2d_cp,tsr,'linear');
    p_z = cp_z.*(1/2).*rho.*(2.*r.*delta_z).*U_inf_zeta.^3;
    ptot = sum(p_z);

end