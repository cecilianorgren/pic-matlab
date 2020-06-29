function ynew = resample_average(x,y,xnew)

nx = numel(xnew);
dx = diff(xnew);
xnew_edges = [x(1), xnew(1:end-1)+0.5*dx, x(end)];


for ix = 1:nx 
  i1 = find(x >= xnew_edges(ix),1,'first'); 
  i2 = find(x <= xnew_edges(ix+1),1,'last'); 
  inds = i1:i2;
  
  ynew(ix) = mean(y(i1:i2));  
end
if 1
  plot(x,y,xnew,ynew)
end