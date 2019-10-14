function out = cross_product(ax,ay,az,bx,by,bz,varargin)
% CROSS_PRODUCT Calculates cross product
% out = cross_product(ax,ay,az,bx,by,bz);
% out = 
% 
%  struct with fields:
% 
%        x: [3200×3200 double]
%        y: [3200×3200 double]
%        z: [3200×3200 double]
%
% out = cross_product(ax,ay,az,bx,by,bz,'components');
% out = 
% 
%  struct with fields:
% 
%        x: [3200×3200 double]
%        y: [3200×3200 double]
%        z: [3200×3200 double]
%     x_yz: [3200×3200 double]
%     x_zy: [3200×3200 double]
%     y_xz: [3200×3200 double]
%     y_zx: [3200×3200 double]
%     z_xy: [3200×3200 double]
%     z_yx: [3200×3200 double]

% Default values
doComponents = 0; % do not return components

% Collect additional inputs
args = varargin;
nargs = numel(varargin);
while not(isempty(args))
  switch args{1}
    case {'components','comps','comp'}
      doComponents = args{2};
      l = 2;
  end
  args(1:l) = [];
end
  
% Calculate cross product
out_x_yz = ay.*bz;
out_x_zy = - az.*by;
out_y_zx = az.*bx;
out_y_xz = - ax.*bz;
out_z_xy = ax.*by;
out_z_yx = - ay.*bx;

out_x = ay.*bz - az.*by;
out_y = az.*bx - ax.*bz;
out_z = ax.*by - ay.*bx;

out.x = out_x;
out.y = out_y;
out.z = out_z;

if doComponents
  if 1
    out.x_yz = out_x_yz;
    out.x_zy = out_x_zy;
    out.y_xz = out_y_xz;
    out.y_zx = out_y_zx;
    out.z_xy = out_z_xy;
    out.z_yx = out_z_yx;
  else    
    out.x.yz = out_x_yz;
    out.x.zy = out_x_zy;
    out.y.xz = out_y_xz;
    out.y.zx = out_y_zx;
    out.z.xy = out_z_xy;
    out.z.yx = out_z_yx;
  end
end