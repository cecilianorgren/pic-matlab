function [new_x,new_z,new_data] = interpolate_to_magnetic_field_arc_length(data,x,z,nq,A,nA)
% [] = interpolate_to_magnetic_field_arc_length(data,x,z,nq,A,nA);
%     
%   1. First interpolates to magnetic field lines obtained from contourcs.
%   2. Interpolate again nq equispaced pounts along magnetic field arc length.
%      
%   Treats top left corner of reconnection region.
%   nA is number of contours or contour levels.
%   

  orig_data = data;
  orig_A = A;
  orig_x = x;
  orig_z = z;
  
  [saddle_locations,saddle_values] = saddle(A,'sort');
  xlim = [x(saddle_locations(1,1)) x(end)];
  zlim = [x(saddle_locations(1,2)) z(end)];
  
  % apply lims first, more efficient 
  [x,z,A,data] = apply_lims_first(x,z,A,data,xlim,zlim);
  S = contourcs(x,z,A',nA);
  lengthLines = [S(:).Length];
  diffLengthLines = diff(lengthLines);
  nLines = numel(S);
  
  doPlot = 1;
  if doPlot
    %%
    h = setup_subplots(3,1);
    isub = 1;
    hca = h(isub); isub = isub + 1;
    imagesc(hca,x,z,A')
    hold(hca,'on')
    for iLine = 1:nLines
      [lineX,lineZ] = apply_lims(S(iLine).X,S(iLine).Y,xlim,zlim);  
      if not(mean(diff(lineX))>1) % increasing X
        lineX = lineX(end:-1:1);
        lineZ = lineZ(end:-1:1);
      end
      plot(hca,lineX,lineZ)
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = zlim;
    
    
    hca = h(isub); isub = isub + 1;
    imagesc(hca,x,z,data')
    hold(hca,'on')
    for iLine = 1:nLines
      [lineX,lineZ] = apply_lims(S(iLine).X,S(iLine).Y,xlim,zlim);  
      plot(hca,lineX,lineZ)
    end    
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = zlim;
    
    hlinks = linkprop(h(1:2),{'XLim','YLim'});
  end
  
  new_z = nan(nq,nLines);
  new_x = nan(nq,nLines);
  new_data = nan(nq,nLines);
  new_lev = nan(1,nLines);
  
  [X,Z] = meshgrid(x,z);
  isnan_line = [];
  isLong = 1;
  for iLine = 1:nLines
    %iLine
    if 0
      lineX = S(iLine).X;
      lineZ = S(iLine).Y;
    else
      [lineX,lineZ] = apply_lims(S(iLine).X,S(iLine).Y,xlim,zlim);  
      levA = S(iLine).Level;
    end
    if isempty(lineX)
      isnan_line = [isnan_line iLine];
      continue
    end
%     if isLong && iLine > 1 && abs(diffLengthLines(iLine-1)) > 0.5*lengthLines(iLine-1) % skip short lines, these are likely islands
%       if isLong == 1
%         isLong = 0;
%       else
%         isLong = 1;
%       end
%     elseif not(isLong) && abs(diffLengthLines(iLine-1)) > 0.5*lengthLines(iLine-1)
%       if isLong == 0
%         isLong = 1;
%       else
%         isLong = 0;
%       end
%     end
    if lengthLines(iLine) > 0.2*max(lengthLines) % simple way, needs to be kept in check manually
      isLong = 1;
    else
      isLong = 0;
    end    
    if not(isLong)
      isnan_line = [isnan_line iLine];
      continue
    end
    Vq = interp2(X,Z,data',lineX,lineZ); % value of data along contour line
    d_arc_x = diff(lineX);
    d_arc_y = diff(lineZ);
    d_arc_distance = sqrt(d_arc_x.^2 + d_arc_y.^2);
    arc_distance = [0 cumsum(d_arc_distance)];
    rel_arc_distance = arc_distance/arc_distance(end);
    new_rel_arc_distance = linspace(rel_arc_distance(1),rel_arc_distance(end),nq);
    
    % resample again from arc_distance to new_rel_arc_distance
    Vq_arc = interp1(rel_arc_distance,Vq,new_rel_arc_distance);
    x_arc = interp1(rel_arc_distance,lineX,new_rel_arc_distance);
    z_arc = interp1(rel_arc_distance,lineZ,new_rel_arc_distance);
    
    if mean(d_arc_x)<0 % if decreasing X, then flip
      x_arc = x_arc(end:-1:1);
      z_arc = z_arc(end:-1:1);
      Vq_arc = Vq_arc(end:-1:1);
    end
    
    new_x(:,iLine) = x_arc;
    new_z(:,iLine) = z_arc;
    new_data(:,iLine) = Vq_arc;   
    new_lev(:,iLine) = levA;
  end
  
  new_x(:,isnan_line) = [];
  new_z(:,isnan_line) = [];
  new_data(:,isnan_line) = [];
  new_lev(:,isnan_line) = [];
  if doPlot
    %%
    hca = h(isub); isub = isub + 1;
    imagesc(hca,new_rel_arc_distance,new_lev,new_data')    
    hold(hca,'on')
    plot(hca,[0 1],saddle_values(1)*[1 1],'k','linewidth',1.5)    
    hold(hca,'off')
    hca.CLim = h(2).CLim;
  end
  
  function [xx,zz] = apply_lims(x,z,xlim,zlim)
    
    xx = x;
    zz = z;
    
    xx(x <= xlim(1)) = NaN;
    zz(x <= xlim(1)) = NaN;
    
    xx(x >= xlim(2)) = NaN;
    zz(x >= xlim(2)) = NaN;
    
    xx(z <= zlim(1)) = NaN;
    zz(z <= zlim(1)) = NaN;
    
    xx(z >= zlim(2)) = NaN;
    zz(z >= zlim(2)) = NaN;
    
    xx(isnan(xx)) = [];
    zz(isnan(zz)) = [];    
    
  end
  function [xx,zz,aa,dd] = apply_lims_first(x,z,A,data,xlim,zlim)
    
    xx = x;
    zz = z;
    aa = A;
    dd = data;
    
    x1 = find(x >= xlim(1),1,'first');
    x2 = find(x <= xlim(2),1,'last');
    z1 = find(z >= zlim(1),1,'first');
    z2 = find(z <= zlim(2),1,'last');
    
    xx = xx(x1:x2);
    zz = zz(z1:z2);    
    aa = aa(x1:x2,z1:z2);
    dd = dd(x1:x2,z1:z2);
    
  end
end