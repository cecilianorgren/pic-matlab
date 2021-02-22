function out = fmax1D(v,n,vd,vt) 
% Get function f(v) for multiple species
%   f = fmax1D(v,n,vd,vt) 

units = irf_units;
nsp = numel(n);
f = @(v,n,vd,vt) n*(1/pi/vt(1)^2)^(1/2)*exp(-(v-vd(1)).^2/vt(1).^2);
ftot = zeros(1,numel(v));
for isp = 1:nsp
  ftot = ftot + f(v,n(isp),vd(isp),vt(isp));
end

out = ftot;