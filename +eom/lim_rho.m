function [value,isterminal,direction] = lim_rho(t,x,X,vx,vy,vz)
% Locate the time when the particle exits the box by left or right and stop
% integration. The limit is adaptive and dependent on the electron energy.

units = irf_units;
energy = 0.5*units.me*(vx.^2+vy.^2+vz.^2);

rho = m*v/q/B;

value = X - abs(x(3))*0.999; % detect z>L/2 and z<L/2 (value = 0)
isterminal = 1;   % Stop the integration
direction = 0;   % Both directions