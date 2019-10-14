function varargout = fun_calc_force_terms(t,x,z,m,q,n,v,p,E,B,varargin)
% FUN_CALC_FORCE_TERMS Calculates force terms

% Default values
doForceDensity = 0; % output is multiplied with density, this reduces large values (incl. inf) in low density regions
doPlot = 0; % plot results, mostly for debugging
doTempdv = 0; % assume only single time input is given
doComponents = 0; % return contributing components of each term, e.g. vxBy

% Collect additional inputs
args = varargin;
nargs = numel(varargin);
while not(isempty(args))
  switch lower(args{1})
    case 'plot'
      doPlot = args{2};
      l = 2;
    case 'density'
      doForceDensity = args{2};
      l = 2;
    case 'comp'
      doComponents = args{2};
      l = 2;
    case {'nmin','minn','nlow','nlim'}
      doNmin = 1;
      nlim = args{2};
      l = 2;
    otherwise
      warning(sprintf('Argument %s not recognized',args{1}));
  end
  args(1:l) = [];
end

% Check if more than one time is given, for calculating dv/dt
if (numel(t) > 1) && (numel(t) == size(v,1))
  
end

%% Calculate force terms
% if time series not given, use given time to load appropriate v in order
% to calculate dv/dt
if size(v,3) > 1 % has several timesteps
  
else
  force_dv_temp.x = v.x*NaN; % not used
  force_dv_temp.y = v.y*NaN; % not used
  force_dv_temp.z = v.z*NaN; % not used
end


dv_conv = convective_derivative(x,z,v,'comp',doComponents); 
vxB = cross_product(v.x,v.y,v.z,B.x,B.y,B.z,'comp',doComponents);
div_p = div_tensor(x,z,p,'comp',doComponents); 
E = struct('x',E.x,'y',E.y,'z',E.z); % sometimes it can have nested structures: E.perp.x

% Perform operation on all fields of structure
% syntax: x = {dv_conv,E,vxB,div_p}
force_dv_conv = structfun(@(x) m*x,  dv_conv ,'UniformOutput',false);
force_E       = structfun(@(x) q*x,  E       ,'UniformOutput',false);
force_vxB     = structfun(@(x) q*x,  vxB     ,'UniformOutput',false);
force_div_p   = structfun(@(x) x./n, div_p   ,'UniformOutput',false);

%% Collect output
varargout{1} = force_dv_temp;
varargout{end+1} = force_dv_conv;
varargout{end+1} = force_E;
varargout{end+1} = force_vxB;
varargout{end+1} = force_div_p;

if doNmin % Remove some moments with very low densities
  indrem = find(n<nlim);  
  for iout = 1:numel(varargout)
    varargout{iout}.x(indrem) = NaN;
    varargout{iout}.y(indrem) = NaN;
    varargout{iout}.z(indrem) = NaN;
  end
end

if doForceDensity % Remove some moments with very low densities    
  for iout = 1:numel(varargout)
    varargout{iout}.x = varargout{iout}.x.*n;
    varargout{iout}.y = varargout{iout}.y.*n;
    varargout{iout}.z = varargout{iout}.z.*n;
  end
end