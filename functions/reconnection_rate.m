function R = reconnection_rate(times,varargin)
% RECONNECTION_RATE Calculates reconnection rate.
%
%   R = reconnection_rate(timesteps/200,'A',A_ts,'E',E_ts(:,:,:,2));
%   
%   
%   plot(timesteps(2:end),R.A,timesteps,R.E)
%   legend({'dA/dt (at X line)','Ey at X line'},'location','best')
%   xlabel('timestep')
%   ylabel('reconnection rate')
%   grid('on')

doA = 1;
doE = 0;

ntimes = numel(times);

nargs = numel(varargin);
args = varargin;

have_option = 1;
while not(isempty(args))
  switch lower(args{1})
    case 'a'
      doA = 1;
      A = args{2};
      l = 2;
    case 'e'
      doE = 1;
      E = args{2};
      l = 2;
  end
  args = args((l+1):end);  
end

if doE
  Ex = zeros(ntimes,1); % electric field at X line
end

Ax = zeros(ntimes,1);
xline_location = zeros(ntimes,2);
for itime = 1:ntimes
  A_tmp = squeeze(A(itime,:,:));
  [saddle_locations,saddle_values] = saddle(A_tmp,'sort');  
  Ax(itime) = saddle_values(1);
  xline_location(itime,1) = saddle_locations(1,1);
  xline_location(itime,2) = saddle_locations(1,2);
  
  if doE
    an_average = 5;
    Ex(itime) = mean(mean(E(itime,saddle_locations(1,1)+an_average*[-1 1],saddle_locations(1,2)+an_average*[-1 1])));
  end
end  

R.Ax = Ax;
R.Xx = xline_location(:,1);
R.Xy = xline_location(:,2);
R.A = diff(Ax)./diff(tocolumn(times));
R.E = Ex;


