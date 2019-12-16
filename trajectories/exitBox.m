function [value, isterminal, direction] = exitBox(T, x)
value      = 200 - abs(x(1)-200); % goes beyond edge of box
isterminal = 1;   % Stop the integration
direction  = 0;