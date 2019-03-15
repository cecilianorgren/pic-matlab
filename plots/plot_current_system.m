% plot_current_system

%% Plot: Parallel vs. perpendicular currents
[X,Z] = ndgrid(x,z); 

xlim = [x(1) x(end)];% + 150*[1 -1];
zlim = [z(1) z(end)]; %zlim = 5*[-1 1];
xlim = [00 70];
zlim = [-6 6];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

doA = 1;
cA = [0 0 0];
doQ = 1;
varstr_Q = 'J';
dataQ = eval(varstr_Q);
nQx = 40;
nQz = 14;
ipxQ = fix(linspace(ipx1,ipx2,nQx));
ipzQ = fix(linspace(ipz1,ipz2,nQz));
ipxQ = ipx1:10:ipx2;
ipzQ = ipz1:10:ipz2;


nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');

clim_jx = 0.1;
clim_jz = 0.1;
clim_by = 0.7;


isub = 1;

% Total current
if 0 % J.par
  varstr = 'J.par';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % J.perp.x
  varstr = 'J.perp.x';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ),'color',[0 0 0]);
    hold(hca,'off')  
  end
end
if 0 % J.perp.y
  varstr = 'J.perp.y';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % J.perp.z
  varstr = 'J.perp.z';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end

% Electron contribution to current
if 0 % -je.par
  varstr = '-je.par';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % -je.perp.x
  varstr = '-je.perp.x';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % -je.perp.y
  varstr = '-je.perp.y';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % -je.perp.z
  varstr = '-je.perp.z';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end

% Electron flux
if 0 % je.par
  varstr = 'je.par';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % je.perp.x
  varstr = 'je.perp.x';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % je.perp.y
  varstr = 'je.perp.y';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % je.perp.z
  varstr = 'je.perp.z';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end

% Ion contribution to current
if 0 % ji.par
  varstr = 'ji.par';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % ji.perp.x
  varstr = 'ji.perp.x';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % ji.perp.y
  varstr = 'ji.perp.y';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 0 % ji.perp.z
  varstr = 'ji.perp.z';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end


if 1 % B.y
  varstr = 'sqrt(ji.x.^2+ji.z.^2)-sqrt(je.x.^2+je.z.^2)';
  %varstr = 'B.z';
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
  hca.CLim = clim_by*[-1 1];
  
  if 1%doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5,'linestyle','--'); 
    hold(hca,'off')  
  end
  if doA
    hold(hca,'on')
    %hcont = contour(hca,x(ipx),z(ipz),ne(ipx,ipz)',[0:0.5:3.5],'color',0*[0.8 0.8 0.8],'linewidth',0.5); 
    % rescale B.y contour to colormap
    byscaling = max(abs(hca.CLim))/max(abs(B.y(:)));
    [ccont,hcont] = contour(hca,x(ipx),z(ipz),B.y(ipx,ipz)'*byscaling,[-0.7:0.05:0.7]*byscaling,'linewidth',1); 
    hold(hca,'off')  
  end
  sQ = 5;
  if 1%doQ
    hold(hca,'on')
    XC = X(ipxQ,ipzQ);
    ZC = Z(ipxQ,ipzQ);
    dataQx = dataQ.x(ipxQ,ipzQ);
    dataQz = dataQ.z(ipxQ,ipzQ);
    %dataQcolor = E.z(ipxQ,ipzQ);%
    dataQcolor = ji.abs(ipxQ,ipzQ)-je.abs(ipxQ,ipzQ); % B.y(ipxQ,ipzQ); % for test
    dataQcolor = sqrt(ji.x(ipxQ,ipzQ).^2+ji.z(ipxQ,ipzQ).^2)-sqrt(je.x(ipxQ,ipzQ).^2+je.z(ipxQ,ipzQ).^2); % B.y(ipxQ,ipzQ); % for test
    nQ = numel(dataQx);
    ind1 = find(dataQcolor>0);
    ind2 = find(dataQcolor<0);
    colors = pic_colors('matlab');
    dataQx1 = dataQx(ind1);
    dataQz1 = dataQz(ind1);
    dataQx2 = dataQx(ind2);
    dataQz2 = dataQz(ind2);
       
    dataQ1 = sqrt(dataQx1.^2+dataQz1.^2); maxQ1 = max(dataQ1(:));
    dataQ2 = sqrt(dataQx2.^2+dataQz2.^2); maxQ2 = max(dataQ2(:));
    maxQ = max([maxQ1,maxQ2]);
    maxQ = 1;
    hquiv1 = quiver(hca,XC(ind1),ZC(ind1),dataQx1*sQ,dataQz1*sQ,0,'color',colors(4,:),'linewidth',3);      
    hquiv2 = quiver(hca,XC(ind2),ZC(ind2),dataQx2*sQ,dataQz2*sQ,0,'color',colors(3,:),'linewidth',3);      
%     for iQ = 1:nQ      
%       if dataQcolor(iQ)>0
%         tmpC = pic_colors('2');
%       else
%         tmpC = pic_colors('4');
%       end
%       %tmpC = 'c';
%       hquiv = quiver(hca,XC(iQ),ZC(iQ),dataQx(iQ),dataQz(iQ),5,'color',tmpC,'linewidth',2);      
%     end
    hold(hca,'off')  
  end
  if 1%doQ
    hold(hca,'on')    
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),J.x(ipxQ,ipzQ)*sQ,J.z(ipxQ,ipzQ)*sQ,0,'color',pic_colors('1'),'linewidth',0.5);
    hold(hca,'off')  
  end
end

h(1).Title.String = sprintf('time = %g (1/wci) = %g (1/wpe)',time,timestep);
hlink = linkprop(h,{'XLim','YLim'});
arrayfun(@(x)eval(sprintf('x.XDir = ''reverse'';'),x),h)
if ncols == 1, compact_panels; end

%% Plot: Main current carriers
for timestep = 5000
  %%
txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/fields-%05.0f.dat',timestep); % michael's perturbation

tic; [x,z,E,B,...
  ni1,ne1,ni2,ne2,...
  vi1,ve1,vi2,ve2,...
  ji1,je1,ji2,je2,...
  pi1,pe1,pi2,pe2,...
  ti1,te1,ti2,te2,...
  dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile); toc
x0 = mean(x); x = x - x0;

% Calculate auxillary quantities
A = vector_potential(x,z,B.x,B.z); % vector potential
[saddle_locations,saddle_values] = saddle(A);
pic_calc_script

varstrs = {'B.y','ve.x','sqrt(ji.x.^2+ji.z.^2)-sqrt(je.x.^2+je.z.^2)'};
nvars = numel(varstrs);

[X,Z] = ndgrid(x,z); 

xlim = [x(1) x(end)];% + 150*[1 -1];
zlim = [z(1) z(end)]; %zlim = 5*[-1 1];
xlim = [00 80];
zlim = [-6 6];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

% Magnetic field lines
doA = 1;
nA = [0:-1:min(A(:))];
cA = [0 0 0];

% Quivers
doQ = 1;
sQ = 5;
varstr_Q = 'J';
dataQ = eval(varstr_Q);
nQx = 40;
nQz = 14;
ipxQ = fix(linspace(ipx1,ipx2,nQx));
ipzQ = fix(linspace(ipz1,ipz2,nQz));
ipxQ = ipx1:10:ipx2;
ipzQ = ipz1:10:ipz2;
%dataQcolor = E.z(ipxQ,ipzQ);%
dataQcolor = ji.abs(ipxQ,ipzQ)-je.abs(ipxQ,ipzQ); % B.y(ipxQ,ipzQ); % for test
dataQcolor = sqrt(ji.x(ipxQ,ipzQ).^2+ji.z(ipxQ,ipzQ).^2)-sqrt(je.x(ipxQ,ipzQ).^2+je.z(ipxQ,ipzQ).^2); % B.y(ipxQ,ipzQ); % for test
indQ{1} = find(dataQcolor>0);
indQ{2} = find(dataQcolor<0);

% Divide current carrying species up into 4
dataQuiverColorMat(:,:,1) = sqrt(ji1.x(ipxQ,ipzQ).^2+ji1.z(ipxQ,ipzQ).^2);
dataQuiverColorMat(:,:,2) = sqrt(je1.x(ipxQ,ipzQ).^2+je1.z(ipxQ,ipzQ).^2);
dataQuiverColorMat(:,:,3) = sqrt(ji2.x(ipxQ,ipzQ).^2+ji2.z(ipxQ,ipzQ).^2);
dataQuiverColorMat(:,:,4) = sqrt(je1.x(ipxQ,ipzQ).^2+je2.z(ipxQ,ipzQ).^2);

colors = pic_colors('matlab');    
colors = colors([4 3 1 2 5 6],:);
    
% Plot
nrows = nvars;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');

clim_jx = 0.1;
clim_jz = 0.1;
clim_by = 0.6;

isub = 1;

for ivar = 1:nvars
  varstr = varstrs{ivar};  
  variable = eval(varstr);
  
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (di)';
  hca.YLabel.String = 'z (di)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.CLim = clim_by*[-1 1];
  
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5,'linestyle','--'); 
    hold(hca,'off')  
  end
  if 1 % By contours
    hold(hca,'on')
    %hcont = contour(hca,x(ipx),z(ipz),ne(ipx,ipz)',[0:0.5:3.5],'color',0*[0.8 0.8 0.8],'linewidth',0.5); 
    % rescale B.y contour to colormap
    byscaling = max(abs(hca.CLim))/max(abs(B.y(:)));
    [ccont,hcont] = contour(hca,x(ipx),z(ipz),B.y(ipx,ipz)'*byscaling,[-0.7:0.05:0.7]*byscaling,'linewidth',1); 
    hold(hca,'off')  
  end
  if doQ % colored by main current carriers
    hold(hca,'on')
    XC = X(ipxQ,ipzQ);
    ZC = Z(ipxQ,ipzQ);
    dataQx = dataQ.x(ipxQ,ipzQ);
    dataQz = dataQ.z(ipxQ,ipzQ);
    nQ = numel(dataQx);
    nindQ = numel(indQ);
    for iindQ = 1:nindQ      
      dataQx_ = dataQx(indQ{iindQ});
      dataQz_ = dataQz(indQ{iindQ});    

      hquiv = quiver(hca,XC(indQ{iindQ}),ZC(indQ{iindQ}),dataQx_*sQ,dataQz_*sQ,0,'color',colors(iindQ,:),'linewidth',3);      
      %hquiv2 = quiver(hca,XC(ind2),ZC(ind2),dataQx2*sQ,dataQz2*sQ,0,'color',colors(3,:),'linewidth',3);
    end      
    hold(hca,'off')  
  end
  if doQ % black
    hold(hca,'on')    
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),J.x(ipxQ,ipzQ)*sQ,J.z(ipxQ,ipzQ)*sQ,0,'color',pic_colors('1'),'linewidth',0.5);
    hold(hca,'off')  
  end
end

h(1).Title.String = sprintf('time = %g (1/wci) = %g (1/wpe)',time,timestep);
hlink = linkprop(h,{'XLim','YLim'});
arrayfun(@(x)eval(sprintf('x.XDir = ''reverse'';'),x),h)

if ncols == 1
  compact_panels; 
  arrayfun(@(x)eval(sprintf('x.XLabel.String = [];'),x),h(1:end-1)); 
end
%cn.print(sprintf('in-plane-current-system_twpe%05.0f',timestep),'path',savedir_root)
end

%% Plot:
[X,Z] = ndgrid(x,z); 

xlim = [x(1) x(end)];% + 150*[1 -1];
zlim = [z(1) z(end)]; %zlim = 5*[-1 1];
xlim = [0 200];
zlim = [0 12];
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


nrows = 3;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');

clim_jx = 0.5;
clim_jz = 0.5;


isub = 1;

if 1 % ji.x
  varstr = 'ji.x';
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
  hca.CLim = clim_jx*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % je.x
  varstr = 'je.x';
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
  hca.CLim = clim_jx*[-1 1];
    
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.x-je.x
  varstr = 'ji.x-je.x';
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
  hca.CLim = clim_jx*[-1 1];
    
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.z
  varstr = 'ji.z';
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
  hca.CLim = clim_jz*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % je.z
  varstr = 'je.z';
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
  hca.CLim = clim_jz*[-1 1];
    
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.z-je.z
  varstr = 'ji.z-je.z';
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
  hca.CLim = clim_jz*[-1 1];
    
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
hlink = linkprop(h,{'XLim','YLim'});
arrayfun(@(x)eval(sprintf('x.XDir = ''reverse'';'),x),h)

%% Plot:
[X,Z] = ndgrid(x,z); 

xlim = [x(1) x(end)];% + 150*[1 -1];
zlim = [z(1) z(end)]; %zlim = 5*[-1 1];
xlim = [0 200];
zlim = [0 12];
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


nrows = 5;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

clim_j = 2;
clim_jJ = 5;

isub = 1;

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
if 1 % ji.x
  varstr = 'ji.x';
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
  hca.CLim = clim_j*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % je.x
  varstr = 'je.x';
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
  hca.CLim = clim_j*[-1 1];
    
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.x
  varstr = 'ji.z';
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
  hca.CLim = clim_j*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % je.x
  varstr = 'je.z';
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
  hca.CLim = clim_j*[-1 1];
    
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.x./J.x
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
  hca.CLim = clim_jJ*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % je.x./J.x
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
  hca.CLim = clim_jJ*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.x./J.x
  varstr = 'ji.x./je.x';
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
  hca.CLim = clim_jJ*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.z./J.z
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
  hca.CLim = clim_jJ*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % je.z./J.z
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
  hca.CLim = clim_jJ*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
if 1 % ji.x./J.x
  varstr = 'ji.z./je.z';
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
  hca.CLim = clim_jJ*[-1 1];
  
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end
hlink = linkprop(h,{'XLim','YLim'});
arrayfun(@(x)eval(sprintf('x.XDir = ''reverse'';'),x),h)