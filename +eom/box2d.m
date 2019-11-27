function [value,isterminal,direction] = box2d(t,x,z,xlim,zlim)
% Locate the time when the particle exits the box by top or bottom and stop
% integration.

isterminal = 1;   % Stop the integration
direction = 0;   % Both directions

if any([x > xlim(2),x < xlim(1),z > zlim(2),z < zlim(1)])
  value = 0;
else 
  value = 1;
end
%value = 1;