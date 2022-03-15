
nPeaks = 10;
spacingPeaks = 2; % decides how far apart the points must be
ftmp = fxyz.f;
v = fxyz.v;
f = fxyz;

for iPeak = 1:nPeaks
  [val,ind] = max(ftmp(:)); % find index and value where ftmp is max

  if val == 0 % skip if phase space density is zero
    continue
  end

  % Change the 'single index' ind to 3 indices, one for each direction
  [ix,iy,iz] = ind2sub(size(ftmp),ind);
  ix_ = ix+[-spacingPeaks:spacingPeaks]; ix_(ix_<0) = [];
  iy_ = iy+[-spacingPeaks:spacingPeaks]; iy_(iy_<0) = [];
  iz_ = iz+[-spacingPeaks:spacingPeaks]; iz_(iz_<0) = [];

  ftmp(ix_,iy_,iz_) = NaN; % remove point and surrounding point so that they are not chosen next time

  % collect results into matrix
  fpeaks(iPeak).vx = f.v(ix);
  fpeaks(iPeak).vy = f.v(iy);
  fpeaks(iPeak).vz = f.v(iz);
  fpeaks(iPeak).x = mean(f.x);
  fpeaks(iPeak).y = 0;
  fpeaks(iPeak).z = mean(f.z);
  fpeaks(iPeak).f = val;

  % just for book keeping
  fpeaks(iPeak).nPeaks = nPeaks;
  fpeaks(iPeak).spacingPeaks = spacingPeaks;           
  fpeaks(iPeak).iPeak = iPeak;
end