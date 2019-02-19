A = vector_potential(x,z,B.x,B.z); % vector potential
KB = magnetic_field_curvature(x,z,B.x,B.y,B.z); % magnetic curvature
b.x = B.x./B.abs;
b.y = B.y./B.abs;
b.z = B.z./B.abs;
% Does this makes no sense, since we deal with the unit vector of B, not the
% absolute value?
% blim = 0.1;
% KB.x(B.abs<blim) = NaN;
% KB.y(B.abs<blim) = NaN;
% KB.z(B.abs<blim) = NaN;
% KB.abs(B.abs<blim) = NaN;

%% Plot
xlim = torow([150 x([round(end/2)])]);
zlim = torow(z([1 end])); zlim = [-5 5];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:2:ipx2;
ipz = ipz1:2:ipz2;

[X,Z] = ndgrid(x,z);

varstr_Q = {'B','KB','KB','KB'};

cQ = [0 0 0];
scaleQ = 1;
ipxQ = ipx1:15:ipx2;
ipzQ = ipz1:15:ipz2; 

cA = 0*[0.8 0.8 0.8];
nA = [0:-0.5:min(A(:))];
ipxA = ipx1:5:ipx2;
ipzA = ipz1:5:ipz2;

doA = 1;
doQ = 1;

% Initialize figure
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end
isub = 1;

if 1 % A
  hca = h(isub); isub = isub + 1;
  varstr = 'A';
  variable = eval(varstr);
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';  
  hcb = colorbar('peer',hca);       
  if doA
    hold(hca,'on')
    hA = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    dataQ = eval(varstr_Q{isub-1});
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ),scaleQ,'color',cQ);
    legend(hquiv,varstr_Q{isub-1})
    hold(hca,'off')  
  end  
  %hca.CLim = [min(abs(himag.CData(:))) max(abs(himag.CData(:)))];
end
if 1 % B.abs
  hca = h(isub); isub = isub + 1;
  varstr = 'B.abs';
  variable = eval(varstr);
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';  
  if 0%any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
    hca.CLim = max(abs(himag.CData(:)))*[-1 0];  
  end
  hcb = colorbar('peer',hca);       
  if doA
    hold(hca,'on')
    hA = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    dataQ = eval(varstr_Q{isub-1});
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ),scaleQ,'color',cQ);
    legend(hquiv,varstr_Q{isub-1})
    hold(hca,'off')  
  end   
  hca.CLim = [min((himag.CData(:))) max((himag.CData(:)))];
end
if 1 % B.abs
  hca = h(isub); isub = isub + 1;
  %imagesc(hca,x(ipx),z(ipz),1./KB.abs(ipx,ipz)')
  varstr = 'KB.abs';
  variable = eval(varstr);
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';  
  if 0%any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
    hca.CLim = max(abs(himag.CData(:)))*[-1 0];  
  end
  hcb = colorbar('peer',hca);   
  if doA
    hold(hca,'on')
    hA = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    dataQ = eval(varstr_Q{isub-1});
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ),scaleQ,'color',cQ);
    legend(hquiv,varstr_Q{isub-1})
    hold(hca,'off')  
  end  
  hca.CLim = [min((himag.CData(:))) max((himag.CData(:)))];
end
if 1 % B.abs
  hca = h(isub); isub = isub + 1;
  %imagesc(hca,x(ipx),z(ipz),1./KB.abs(ipx,ipz)')
  varstr = 'KB.z';
  variable = eval(varstr);
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';  
  if 0%any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
    hca.CLim = max(abs(himag.CData(:)))*[-1 0];  
  end
  hcb = colorbar('peer',hca);   
  if doA
    hold(hca,'on')
    hA = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    dataQ = eval(varstr_Q{isub-1});
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ),scaleQ,'color',cQ);
    legend(hquiv,varstr_Q{isub-1})
    hold(hca,'off')  
  end  
  hca.CLim = [min((himag.CData(:))) max((himag.CData(:)))];
end
c_eval('colormap(h(?),pic_colors(''candy''));',1:numel(h))
c_eval('axis(h(?),''equal'');',1:4)
c_eval('h(?).XLim = xlim;',1:4)
c_eval('h(?).YLim = zlim;',1:4)
c_eval('h(?).Position(3) = h(1).Position(3);',1:4)
c_eval('h(?).Position(4) = h(1).Position(4);',1:4)