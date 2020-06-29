function out = boxcar_downsampling(t,x,width,step)

n = numel(t);
start = 1:step:(n-step);
stop = start + width;
stop(stop > n) = n;
for iwin = 1:numel(start)
  newt(iwin) = mean(t(start(iwin):stop(iwin)));
  newx(iwin) = mean(x(start(iwin):stop(iwin)));
end
if 1
  plot(t,x,newt,newx)
end