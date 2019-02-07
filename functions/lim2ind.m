function out = lim2ind(varargin)
if nargin == 2  
  val = varargin{1};
  lim = varargin{2};
  ifirst = find(val>lim(1),1,'first');
  ilast = find(val<lim(2),1,'last');
  out = ifirst:1:ilast;
elseif nargin == 4
  val1 = varargin{1};
  lim1 = varargin{2};
  val2 = varargin{3};
  lim2 = varargin{4};
  
%   ifirst1 = find(val1>lim1(1),1,'first');
%   ilast1 = find(val1<lim1(2),1,'last');
%   ifirst2 = find(val2>lim2(1),1,'first');
%   ilast2 = find(val2<lim2(2),1,'last');
%   ind1 = ifirst1:ilast1;
%   ind2 = ifirst2:ilast2;
  
  [VAL1,VAL2] = meshgrid(val1,val2);
  iFirst1 = find(VAL1>lim1(1));
  iLast1 = find(VAL1<lim1(2));
  i1 = intersect(iFirst1,iLast1);
  iFirst2 = find(VAL2>lim2(1));
  iLast2 = find(VAL2<lim2(2));
  i2 = intersect(iFirst2,iLast2);
  out = intersect(i1,i2);  
end