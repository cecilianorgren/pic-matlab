function out = f_maxwellian(n,vd,vt,dim) 
% Get function f(v) for multiple species
%   f = f_maxwellian(n,vd,vt,dim) 
%   where f is a function f(v1), f(v1,v2), or f(v1,v2,v3), as decided by
%   the argument 'dim' 1,2, or 3, the dimension of vd and vt must match dim

units = irf_units;
nsp = numel(n);

switch dim
  case 1
    f = @(vx,n,vd,vt) n*(1/pi/vt(1)^2)^(1/2)*exp(-(vx-vd(1)).^2/vt(1).^2);
    ftot = @(vx) 0;
    for isp = 1:nsp
      ftot = @(vx) ftot(vx) + f(vx,n(isp),vd(isp),vt(isp));
    end
  case 2
    f = @(vx,vy,n,vd,vt) n*(1/pi/vt(1)^2)^(2/2)*exp(-(vx-vd(1)).^2/vt(1).^2-(vy-vd(2)).^2/vt(2).^2);  
    ftot = @(vx,vy) 0;
    for isp = 1:nsp
      ftot = @(vx,vy) ftot(vx,vy) + f(vx,vy,n(isp),vd(isp),vt(isp));
    end    
  case 3
    vttot = sqrt(sum(vt.^2));
    f = @(vx,vy,vz,n,vd,vt) n*(1/pi/vttot^2)^(3/2)*exp(-(vx-vd(1)).^2/vt(1).^2-(vy-vd(2)).^2/vt(2).^2-(vz-vd(3)).^2/vt(3).^2);    
    ftot = @(vx,vy,vz) 0;
    for isp = 1:nsp
      ftot = @(vx,vy,vz) ftot(vx,vy,vz) + f(vx,vy,vz,n(isp),vd(isp),vt(isp));
    end    
end

out = ftot;