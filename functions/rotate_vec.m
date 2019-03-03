function vargout = rotate_vec(varargin)

switch nargin 
  case 2
    old_x = varargin{1}.x;
    old_y = varargin{1}.y;
    old_z = varargin{1}.z;        
    rx = varargin{2}(1,:);
    ry = varargin{2}(2,:);
    rz = varargin{2}(3,:);
  case 4
    old_x = varargin{1}.x;
    old_y = varargin{1}.y;
    old_z = varargin{1}.z;
    rx = varargin{2};
    ry = varargin{3};
    rz = varargin{4};   
  case 6
    old_x = varargin{1};
    old_y = varargin{2};
    old_z = varargin{3};
    rx = varargin{4};
    ry = varargin{5};
    rz = varargin{6};      
  otherwise
    error('Input not recognized.')    
end

new_x = old_x.*rx.x + old_y.*rx.y + old_z.*rx.z;
new_y = old_x.*ry.x + old_y.*ry.y + old_z.*ry.z;
new_z = old_x.*rz.x + old_y.*rz.y + old_z.*rz.z;

switch nargout
  case 1
    new_xyz.x = new_x;
    new_xyz.y = new_y;
    new_xyz.z = new_z;
    vargout(1) = new_xyz;
  case 3
    vargout(1) = new_x;
    vargout(2) = new_y;
    vargout(3) = new_z;
end
