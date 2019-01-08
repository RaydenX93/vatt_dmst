function ret = virtual(npoints, omega, c, x_from, r, alpha, U0, U1, sinTh, cosTh)
% Evaluates the virtual camber correction to be added to the geometric angle of attack
% ----------------------------------------------
% INPUT
% ----------------------------------------------
% npoints = number of virtual points along chord line on which the virtual
% correction will be evaluated
% omega = rotational velocity in rad/s
% c = airfoil chord
% x_from = x_from_nose declared in main routine
% r = cell radius
% alpha = geometric angle of attack
% U0 = local absolute x-velocity
% U1 = local absolute y-velocity
% sinTh = local azimuthal angle sine
% cosTh = local azimuthal angle cosine
% ----------------------------------------------
% OUTPUT
% ----------------------------------------------
% ret = virtual camber correction to be added to the geometric
% angle of attack

% local reference system:
% origin in airfoil pivoting point
% x positive towards trailing edge
% y positive along outward radial direction

a = 0;
V1 = -U0.*sinTh+U1.*cosTh;
V2 = U0.*cosTh+U1.*sinTh;

for i=0:npoints-1
   pos_x = (i.*c)./(npoints-1)-x_from;
   r_prime = sqrt(r.^2+pos_x.^2);
   delta=atan(pos_x./r);
   
   l1 = V1 - omega.*r_prime.*sin(delta); % perpendicular to chord %
   l2 = V2 + omega.*r_prime.*cos(delta); % parallel to chord %
   
   alpha_prime = radtodeg(atan( - l1./ l2)); % minus due to our convention of AoAs %
   local_virt = alpha_prime - alpha;
   
   a = a + local_virt;
end

ret = a./npoints;
end