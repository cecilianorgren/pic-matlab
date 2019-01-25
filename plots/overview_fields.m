%% Load data
timestep = 04800;
txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons/data/fields-%05.0f.dat',timestep); % michael's perturbation

%timestep = 055;
%txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons_2/data/fields-%05.0f.dat',timestep); % try-out with larger perturbation


%tic; [time,r,e,b,n1,ne,ve,vi,je,ji,pe,pi,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q] = read_fields_ieie(txtfile); toc
tic; [x,z,E,B,...
  ni1,ne1,ni2,ne2,...
  vi1,ve1,vi2,ve2,...
  ji1,je1,ji2,je2,...
  pi1,pe1,pi2,pe2,...
  ti1,te1,ti2,te2,...
  dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile); toc

% Calculate auxillary quantities
tic; A = vector_potential(x,z,B.x,B.z); toc % vector potential
pb = B.abs.^2/2; % magnetic pressure
bcurv = magnetic_field_curvature(x,z,B.x,B.y,B.z); % magnetic curvature
c_eval('ve?xB = cross_product(ve?.x,ve?.y,ve?.z,B.x,B.y,B.z);',1:2) % electron motional electric field
c_eval('vi?xB = cross_product(vi?.x,vi?.y,vi?.z,B.x,B.y,B.z);',1:2) % ion motional electric field
ExB = cross_product(E.x,E.y,E.z,B.x,B.y,B.z); % Poynting flux
c_eval('E_ve?xB.x = E.x + ve?xB.x; E_ve?xB.y = E.y + ve?xB.y; E_ve?xB.z = E.z + ve?xB.z;',1:2) % electron motional electric field
c_eval('E_vi?xB.x = E.x + vi?xB.x; E_vi?xB.y = E.y + vi?xB.y; E_vi?xB.z = E.z + vi?xB.z;',1:2) % electron motional electric field
c_eval('je?E = je?.x.*E.x + je?.y.*E.y + je?.y.*E.z;',1:2)
c_eval('ji?E = ji?.x.*E.x + ji?.y.*E.y + ji?.y.*E.z;',1:2)
UB = 0.5*B.abs.^2;
c_eval('Uek? = mass(2)*0.5*ne?.*(ve?.x.^2 + ve?.y.^2 + ve?.z.^2);',1:2)
c_eval('Uik? = mass(1)*0.5*ni?.*(vi?.x.^2 + vi?.y.^2 + vi?.z.^2);',1:2)
c_eval('Uet? = pe?.scalar;',1:2)
c_eval('Uit? = pi?.scalar;',1:2)
Uktot = Uik1 + Uik2 + Uek1 + Uek2;
Uttot = Uit1 + Uit2 + Uet1 + Uet2;
jtot.x = ji1.x + ji2.x - je1.x - je2.x;
jtot.y = ji1.y + ji2.y - je1.y - je2.y;
jtot.z = ji1.z + ji2.z - je1.z - je2.z;

%% Plot 1, 4 species plasma properties, 1 species per column
% Initialize figure
nrows = 5;
ncols = 4;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end

% Panels
doA = 0;
cA = [0.8 0.8 0.8];
nA = 20;
nA = [0:-2:min(A(:))];
ipx = 1:2:nnx;
ipz = 1:2:nnz;
isub = 1;
if 0 % A
  hca = h(isub); isub = isub + 1;
  varstr = 'A';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr; 
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  hcb.YLim = hca.CLim(2)*[-1 0];
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % babs
  hca = h(isub); isub = isub + 1;
  varstr = 'B.abs';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr; 
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  hcb.YLim = hca.CLim(2)*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % bx
  hca = h(isub); isub = isub + 1;
  varstr = 'B.x';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % by
  hca = h(isub); isub = isub + 1;
  varstr = 'B.y';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % bz
  hca = h(isub); isub = isub + 1;
  varstr = 'B.z';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % epar
  hca = h(isub); isub = isub + 1;
  varstr = 'E.par';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);  
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ex
  hca = h(isub); isub = isub + 1;
  varstr = 'E.x';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ey
  hca = h(isub); isub = isub + 1;
  varstr = 'E.y';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ez
  hca = h(isub); isub = isub + 1;
  varstr = 'E.z';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ne1
  hca = h(isub); isub = isub + 1;
  varstr = 'ne1';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ne2
  hca = h(isub); isub = isub + 1;
  varstr = 'ne2';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ni1
  hca = h(isub); isub = isub + 1;
  varstr = 'ni1';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ni2
  hca = h(isub); isub = isub + 1;
  varstr = 'ni2';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vepar1
  hca = h(isub); isub = isub + 1;  
  varstr = 've1.par';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;  
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);  
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vepar2
  hca = h(isub); isub = isub + 1;
  varstr = 've2.par';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);  
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vipar1
  hca = h(isub); isub = isub + 1;
  varstr = 'vi1.par';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);  
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vipar2
  hca = h(isub); isub = isub + 1;
  varstr = 'vi2.par';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);  
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vex1
  hca = h(isub); isub = isub + 1;
  varstr = 've1.x';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vex2
  hca = h(isub); isub = isub + 1;
  varstr = 've2.x';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vix1
  hca = h(isub); isub = isub + 1;
  varstr = 'vi1.x';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vix2
  hca = h(isub); isub = isub + 1;
  varstr = 'vi2.x';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vey1
  hca = h(isub); isub = isub + 1;
  varstr = 've1.y';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vey2
  hca = h(isub); isub = isub + 1;
  varstr = 've2.y';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % viy1
  hca = h(isub); isub = isub + 1;
  varstr = 'vi1.y';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % viy2
  hca = h(isub); isub = isub + 1;
  varstr = 'vi2.y';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vez1
  hca = h(isub); isub = isub + 1;
  varstr = 've1.z';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vez2
  hca = h(isub); isub = isub + 1;
  varstr = 've2.z';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % viz1
  hca = h(isub); isub = isub + 1;
  varstr = 'vi1.z';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % viz2
  hca = h(isub); isub = isub + 1;
  varstr = 'vi2.z';
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = varstr;
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pxx e1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe1.xx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e1xx}';
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pxx e2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe2.xx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e2xx}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pxx i1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi1.xx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i1xx}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pxx i2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi2.xx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i2xx}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pyy e1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe1.yy');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e1yy}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pyy e2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe2.yy');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e2yy}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pyy i1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi2.yy');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i1yy}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pyy i2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi1.yy');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i2yy}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pzz e1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe1.zz');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e1zz}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pzz e2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe2.zz');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e2zz}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pzz i1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi1.zz');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i1zz}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pzz i2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi2.zz');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{izz2}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pscalar e1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe1.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e1scalar}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pscalar e2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pe2.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{e2scalar}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pscalar i1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi1.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i1scalar}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % pscalar i2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,pi2.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'p_{i2scalar}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb.YLim = [0 hca.CLim(2)];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ve1xB_x
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ve1xB.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{e1}xB)_x';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ve2xB_x
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ve2xB.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{e2}xB)_x';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vi1xB_x
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,vi1xB.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{i1}xB)_x';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vi2xB_x
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,vi2xB.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{i2}xB)_x';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ve1xB_y
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ve1xB.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{e1}xB)_y';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ve2xB_y
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ve2xB.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{e2}xB)_y';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vi1xB_y
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,vi1xB.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{i1}xB)_y';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vi2xB_y
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,vi2xB.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{i2}xB)_y';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ve1xB_z
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ve1xB.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{e1}xB)_z';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ve2xB_z
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ve2xB.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{e2}xB)_z';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vi1xB_z
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,vi1xB.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{i1}xB)_z';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % vi2xB_z
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,vi2xB.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '(v_{i2}xB)_z';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % je1E
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,je1E');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{e1}\cdot E';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % je2E
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,je2E');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'je2E';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ji1E
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ji1E');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'ji1E';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % ji2E
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,ji2E');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'ji1E';
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uek1
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uek1');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{ek1}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uek2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uek2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{ek2}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uik1
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uik1');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{ik1}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uik2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uik2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{ik2}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uet1
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uet1');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{et1}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uet1
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uet2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{et2}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uit1
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uit1');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{it1}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % Uit2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uit2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{it2}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % Uktot
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uktot');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{ktot}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % Uktot
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,Uttot');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{ttot}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % UB
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,UB');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{B}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 0 % UB + Uk + Ut
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,UB' + Uttot' + Uktot');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'U_{B} + U_{ttot} + U_{ktot}';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hcb.YLim(1) = 0;
  hca.CLim = hca.CLim(2)*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end


for ipanel = 1:npanels
  h(ipanel).YDir = 'normal';
  h(ipanel).YLim = [-10 10];
end

%% Plot, define variable in cell array
% Define what variables to plot
varstrs = {'ve1.x','ve2.x','ve1.z','ve2.z','ve1.par','ve2.par','ExB.x','ExB.z','-ve1xB.x','-ve2xB.x','-ve1xB.z','-ve2xB.z','E.x','E.z'};
nvars = numel(varstrs);

% Initialize figure
npanels = nvars;
nrows = 7;
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end

% Panels
doA = 0;
cA = [0.8 0.8 0.8];
nA = 20;
nA = [0:-2:min(A(:))];
ipx = 1:2:nnx;
ipz = 1:2:nnz;
isub = 1;
for ivar = 1:nvars  
  hca = h(isub); isub = isub + 1;
  varstr = varstrs{ivar};
  variable = eval(varstr);  
  himag = imagesc(hca,x,z,variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  hcb = colorbar('peer',hca);
  %hcb.YLim = hca.CLim(2)*[-1 1];
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end


for ipanel = 1:npanels
  h(ipanel).YDir = 'normal';
  h(ipanel).XLim = h(ipanel).XLim + [140 -140];
  h(ipanel).YLim = [-5 5];
end