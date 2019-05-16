function varargout = dAdt(t,A,varargin)

if not(isempty(varargin)) && strcmp(varargin{1},'xline'); doXline = 1; else doXline = 0; end

varargin = varargin(2:end);
if not(isempty(varargin)) && strcmp(varargin{1},'xline'); B = varargin{2}; end

matsize = size(A);
dAdt = zeros(matsize); % partial time derivative

dt = t(2)-t(1);

dAdt(1,:,:) = (A(2,:,:)-A(1,:,:))/dt;
dAdt(2:end-1,:,:) = (A(3:end,:,:)-A(1:end-2,:,:))/(2*dt);
dAdt(end,:,:) = (A(end,:,:)-A(end-1,:,:))/dt;

out = dAdt;


%   if doXline
%     dAdt_xline = zeros(numel(t),1);
%     for it = 1:numel(t)
%       [saddle_locations,saddle_values] = saddle(squeeze(A(it,:,:)),'sort');
%       dAdt_xline(it) = dAdt(it,saddle_locations(1,1),saddle_locations(1,2));
%     end
%     out = dAdt_xline;
%   end

  if doXline
    dAdt_xline = zeros(numel(t),1);
    A_saddle = zeros(numel(t),1);
    for it = 1:numel(t)
      [saddle_locations,saddle_values] = saddle(squeeze(A(it,:,:)),'sort');
      A_saddle(it) = saddle_values(1);    
    end
    dAdt_saddle_tmp = diff(A_saddle)/dt;
    dAdt_saddle = interp1(t(1:end-1)+0.5*dt,dAdt_saddle_tmp,t);
    out = dAdt_saddle;
  end

if nargout == 1
  varargout{1} = out;
elseif nargout == 2
  varargout{1} = out;
  varargout{2} = A_saddle;
end
end % end of main function