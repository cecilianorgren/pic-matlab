function out = ind_from_lim(val,lim,varargin)
% IND_FROM_LIM find indices inside from x and z values
%   inds = ind_from_lim(val,lim,valargin)
%   val - xe,ze, xi, or zi (the lims has to be in the same units)
%   lim - xlim, or zlim - one or two values each
%                if you use single values, you might want to supply an
%                option
%   Options:
%     'bounding',nBounding - find a given number of values bounding the given value
%     'closest',nClosest - find the given number of closest values
%
%   Examples: 
%     % Get indices within a range of values
%     ind = ind_from_lim(xi,[18 20]);
%     xi(ind)
%   
%     % Get n number of bounding values
%     % (this should give the same as example above)
%     ind = ind_from_lim(xi,20,'bounding',3);
%     xi(ind)

%     % Get n number of indices closest to a given value  
%     ind = ind_from_lim(xi,20,'closest',5);
%     xi(ind)
%
%     % Plot subpart of data
%     xind = ind_from_lim(xi,[190 210]);
%     zind = ind_from_lim(zi,[-4 4]);
%     pcolor(xi(xind),zi(zind),Bx(xind,zind)')
%     
%     % Pick out a single value
%     xind = ind_from_lim(xi,204);
%     zind = ind_from_lim(zi,2);
%     Bx(xind,zind)
%     
%     % Pick out a single value from a small area around given pount
%     xind = ind_from_lim(xi,204,'closest',5);
%     zind = ind_from_lim(zi,2,'closest',5);
%     mean(mean(Bx(xind,zind)))

% Defaults
doExact = 0; % find all exact matches, can be any number
doBounding = 0; % find a given number of values bounding the given value
nBounding = 0;
doClosest = 0;
nClosest = 1; % only the closest index, can be one or many, for example to reduce cadence

if numel(val) == 1
  doClosest = 1;        
end

% Check input
have_options = 0;
nargs = numel(varargin);      
if nargs > 0, have_options = 1; args = varargin(:); end

while have_options
  l = 1;
  switch(lower(args{1}))
    case 'closest'
      l = 2;
      doClosest = 1;
      nClosest = args{2};            
      args = args(l+1:end);    
    case 'bounding'
      l = 2;
      doBounding = 1;
      nBounding = args{2};
      args = args(l+1:end);
    case 'exact'
      l = 1;
      doExact = 1;
      args = args(l+1:end);
    otherwise
      warning(sprintf('Input ''%s'' not recognized.',args{1}))
      args = args(l+1:end);
  end        
  if isempty(args), break, end    
end

% Find indices
if doBounding
  i0 = find(abs(val-lim(1)) == min(abs(val-lim(1))));
  i1 = i0 - nBounding;
  i2 = i0 + nBounding;
  % Check so that indices are not outside range
  if i1 < 1
    i1 = 1; 
  end
  if i2 > numel(val) 
    i2 = numel(val); 
  end
  inds = i1:i2;
elseif doClosest        
  ii = abs(val-lim(1));
  [is, index] = sort(abs(val-lim(1)));
  inds = sort(index(1:nClosest)); % why did i just not take the first closest index?     
elseif doExact
  [~,inds,~] = intersect(val,lim);
else        
  i1 = find(val >= lim(1),1,'first'); 
  i2 = find(val <= lim(2),1,'last'); 
  inds = i1:i2;
end

out = inds;
