function [value,isterminal,direction] = lim(t,x,X)
% Locate the time when the particle exits the box by top or bottom and stop
% integration.

value = X - abs(x(3))*0.999; % detect z>L/2 and z<L/2 (value = 0)
isterminal = 1;   % Stop the integration
direction = 0;   % Both directions