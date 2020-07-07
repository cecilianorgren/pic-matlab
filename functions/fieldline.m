function varargout = fieldline(x0,z0,x,z,Bx,By,Bz,dx,dy,dz,nsteps)
% FIELDLINE Integrates field in direction of flow vector
%   Example: 
%   x0 = 150;
%   z0 = 9.6;
%   dx = 0.1;
%   dy = 0.1;
%   dz = 0.1;
%   nsteps = 700; % if nsteps is not an even number, integration will stop 
%                 % when the fieldline arclength is above lwngth defined by
%                 % nsteps (it's a bit ugly)
%   nsteps = 100.1;
%   [linearclength,linex,liney,linez,linebx,lineby,linebz] = fieldline(x0,z0,x,z,Bx,By,Bz,dx,dy,dz,nsteps);

  % If nsteps is an uneven function (checked below), process is stopped
  % when line arclength is above nsteps
  doArclengthLimit = 0;
  % Plot progress (for debugging)
  doPlot = 0; 
  % Make field unit vectors
  babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
  bx = Bx./babs;
  by = By./babs;
  bz = Bz./babs;
  
  % Set default values if not given in input
  if isempty(dx); dx = 0.2; end
  if isempty(dy); dy = 0.2; end
  if isempty(dz); dz = 0.2; end
  if isempty(nsteps); nsteps = 100; end
  if not(mod(nsteps,1)==0); doArclengthLimit = 1; end
  
  % Define starting point
  arcline = 0;
  xline = x0;  
  yline = 0;
  zline = z0;  
  % Interploate the magnetic field to starting point
  bxline = interpfield(x,z,bx,xline(end),zline(end));
  byline = interpfield(x,z,by,xline(end),zline(end));
  bzline = interpfield(x,z,bz,xline(end),zline(end));
  
  if doPlot % Plot starting point
    plot(xline(end),zline(end),'.')
    drawnow;
    hold(gca,'on')
  end

  istep = 0;
  while 1
    istep = istep + 1;
    % Calculate step, magnetic field direction times the predefined stepsize
    xstep = bxline(istep)*dx;
    ystep = byline(istep)*dy;
    zstep = bzline(istep)*dz;    
    % Advance line
    xline(istep+1) = xline(istep) + xstep;
    yline(istep+1) = yline(istep) + ystep;
    zline(istep+1) = zline(istep) + zstep;
    arcline(istep+1) = arcline(istep) + sqrt(xstep^2+ystep^2+zstep^2);
    
    % Break if we end up outside of box
    if or(xline(end)<x(1),xline(end)>x(end))
      xline(end) = [];
      yline(end) = [];
      zline(end) = [];
      break;
    end
    % Get field at new point    
    bxline(istep+1) = interpfield(x,z,bx,xline(end),zline(end));
    byline(istep+1) = interpfield(x,z,by,xline(end),zline(end));
    bzline(istep+1) = interpfield(x,z,bz,xline(end),zline(end));
      
    if doPlot
      plot(xline(end),zline(end),'.')
      drawnow;
    end
    
    if doArclengthLimit 
      if arcline(end)>nsteps
        break;
      end
    elseif istep > (nsteps-1)
      break;
    end
    
  end
  if doPlot
    hold(gca,'off')
  end
  
  
  varargout{1} = arcline;
  varargout{2} = xline;
  varargout{3} = yline;
  varargout{4} = zline;
  varargout{5} = bxline;
  varargout{6} = byline;
  varargout{7} = bzline;
end