   
  function [value, isterminal, direction] = myEventBoxEdge(t, x_vect,boxedge)
    % integration is terminated when value changes sign
    % for this setup, value is initially negative
    value      = (x_vect(1)-boxedge(1))*(x_vect(1)-boxedge(2));        
    isterminal = 1;   % Stop the integration
    direction  = 0;
  end