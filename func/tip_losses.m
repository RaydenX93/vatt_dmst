function ret = tip_losses(blades, mu_z, senAlpha)
% Evaluates the tip loss factor along the blade
% ----------------------------------------------
% INPUT
% ----------------------------------------------
% blades = number of blades of the turbine
% mu_z = non-dimensional vertical position
% senAlpha = local angle of attack
% ----------------------------------------------
% OUTPUT
% ----------------------------------------------
% ret = non-dimensional tip loss factor

f1 = abs(mu_z);

% Prandtl-Glauert formula %
if abs(senAlpha) < eps
    senAlpha = 0.001;
end

g1 = exp((1*blades*(f1-1))/(2*f1*abs(senAlpha)));
ret = (2/pi)*acos(min([1 g1]));

end