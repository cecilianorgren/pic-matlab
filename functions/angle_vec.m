function out = angle_vec(varargin)
% Calculates angle between vectors
switch nargin
  case 2
    vec1 = varargin{1};
    vec2 = varargin{2};
    x1 = vec1.x;
    y1 = vec1.y;
    z1 = vec1.z;
    x2 = vec2.x;
    y2 = vec2.y;
    z2 = vec2.z;
    
    abs1 = sqrt(x1.^2 + y1.^2 + z1.^2);
    abs2 = sqrt(x2.^2 + y2.^2 + z2.^2);
    
    dot_prod = (x1.*x2 + y1.*y2 + z1.*z2)./abs1./abs2;

    angle = acosd(dot_prod);
    out = angle;
  case 4
    
  case 6
    
  otherwise
    error('Wrong number of inputs.')
end  
