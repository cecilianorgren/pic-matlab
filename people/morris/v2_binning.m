function [f,v2mid,vol] = v2_binning(v,F,v2edges)

dv_orig = v(2)-v(1);
[VX,VY,VZ] = ndgrid(v,v,v);
V2 = VX.^2 + VY.^2 + VZ.^2;

V2 = V2(:);
F = F(:);

vedges = sqrt(v2edges); % v = sqrt(v^2)
vmid = vedges(2:end)-0.5*(vedges(2:end)-vedges(1:end-1));
dv = diff(vedges);
vol = vmid.^2.*dv;
vol = vol';

[NF edges mid loc] = histcn(V2,v2edges,'AccumData',F*dv_orig.^3);

f = NF./vol;
v2mid = mid{1};

