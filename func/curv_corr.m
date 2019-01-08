function f = curv_corr(tsr, T3)
% Streamtube curvature correction %
pol_c = [-1.81092769088177,2.20489630426131,0.0200817334712570,-0.750931561119896,-0.0341679415246660,0.000121233546060024,0.0820475112021503,0.0160284941069508,8.33599248210428e-05,-2.35310440462014e-06,-0.00195957452056302,-8.76961301816267e-05,6.37387403753171e-07,7.57339447006923e-09,1.14200105541618e-05,-1.20128313460656e-08,-6.35820934586062e-10,-7.00204673074453e-12];

if tsr < 2.1
    tsr = 2.1;
elseif tsr > 3.6
    tsr = 3.6;
end

f = pol_c(1) + pol_c(2)*tsr + pol_c(3)*rad2deg(T3) + pol_c(4)*tsr^2 + pol_c(5)*tsr*rad2deg(T3) + pol_c(6)*rad2deg(T3)^2 + pol_c(7)*tsr^3 + pol_c(8)*tsr^2*rad2deg(T3) + pol_c(9)*tsr*rad2deg(T3)^2 + pol_c(10)*rad2deg(T3)^3 + pol_c(11)*tsr^3*rad2deg(T3) + pol_c(12)*tsr^2*rad2deg(T3)^2 + pol_c(13)*tsr*rad2deg(T3)^3 + pol_c(14)*rad2deg(T3)^4 + pol_c(15)*tsr^3*rad2deg(T3)^2 + pol_c(16)*tsr^2*rad2deg(T3)^3 + pol_c(17)*tsr*rad2deg(T3)^4 + pol_c(18)*rad2deg(T3)^5;
end