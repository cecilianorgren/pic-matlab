%phi = scalar_potential(x,z,E.x,E.z); % vector potential
phi = phi_orig(1:end-1,1:end-1);

[X,Z] = ndgrid(x,z);
  
if 0
  dXx = X(2:end,:)-0.5*diff(X,1);
  dn = ni-ne;
  dphi_x = diff(phi,1);
  phi_Ez = interp2(dXx,dphi_x,X);
end

xlim = [x(1) x(end)];% + 150*[1 -1];
zlim = [z(1) z(end)]; %zlim = 5*[-1 1];
xlim = [-200 0];
zlim = [0 15];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

doQ = 1;
varstr_Q = 'J';
dataQ = eval(varstr_Q);
nQx = 40;
nQz = 10;
ipxQ = fix(linspace(ipx1,ipx2,nQx));
ipzQ = fix(linspace(ipz1,ipz2,nQz));
ipxQ = ipx1:20:ipx2;
ipzQ = ipz1:20:ipz2;


nrows = 2;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

isub = 1;
if 0 % phi
  varstr = 'phi';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % Ez
  varstr = 'E.z';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % Ex
  varstr = 'E.x';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % log2(ji./je).x
  varstr = 'ji.x./J.x';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % log2(ji./je).x
  varstr = 'je.x./J.x';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % log2(ji./je).x
  varstr = 'ji.z./J.z';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % log2(ji./je).x
  varstr = 'je.z./J.z';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % diff(phi,1)
  varstr = 'diff(phi,2)';
  variable = eval(varstr);
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,variable');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % ni-ne
  varstr = 'ni-ne';
  hca = h(isub); isub = isub + 1;
  variable = eval(varstr);
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = input_varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % diff(phi)/dx
  varstr = 'phi';
  hca = h(isub); isub = isub + 1;
  variable = eval(varstr);
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = input_varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
hlink = linkprop(h,{'XLim','YLim'});
