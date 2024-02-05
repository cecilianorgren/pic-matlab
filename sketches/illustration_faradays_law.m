B0 = 0.5;
L = 1;

E0 = 0.2;
Lx = 1;
Lz = 1;
t1 = 1;
t2 = 2.1;

Bx = @(z) B0*tanh(z/L);
Ey = @(x,z,t) E0*exp(-x.^2/(Lx*t)^2-z.^2/(Lz*t)^2);


x = linspace(-4,4,100); dx = x(2) - x(1);
z = linspace(-2,2,40); dz = z(2) - z(1);
[X,Z] = ndgrid(x,z);

BX = Bx(Z);
EY = Ey(X,Z,t1);
EY1 = Ey(X,Z,t1);
EY2 = Ey(X,Z,t2);
BZ = Z*0;
%EY = Ey(X,Z,t2);

rotEY_x = -[diff(EY,1,2)/dz, (EY(:,end)-EY(:,end-1))/dz];
rotEY_z = [diff(EY,1,1)/dz; (EY(end,:)-EY(end-1,:))/dz];
%dEYdt = (EY2-EY1)/(t2-t1);

% For quivers.
xq = linspace(-4,4,20); dxq = xq(2) - xq(1);
zq = linspace(-2,2,15); dzq = zq(2) - zq(1);
[Xq,Zq] = ndgrid(xq,zq);


BXq = Bx(Zq);
EYq = Ey(Xq,Zq,t1);
BZq = Zq*0;
%EY = Ey(X,Z,t2);

rotEYq_x = -[(EYq(:,2)-EYq(:,1))/dzq, (EYq(:,3:end)-EYq(:,1:end-2))/(2*dzq), (EYq(:,end)-EYq(:,end-1))/dzq];
rotEYq_z = [(EYq(2,:)-EYq(1,:))/dxq; (EYq(3:end,:)-EYq(1:end-2,:))/(2*dxq); (EYq(end,:)-EYq(end-1,:))/dxq];

h = setup_subplots(1,1);
isub = 1;
cmap = pic_colors('blue_red'); cmap = flipdim(cmap(1:floor(size(cmap,1)/2),:),1);
colors = pic_colors('matlab');

if 0
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,BX)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'B_x';
end
if 0
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,EY1)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'E_y(t_1)';
end
if 0 % EY2
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,EY2)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_y(t_2)';
end
if 0 % dEYdt
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,dEYdt)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'dE_y/dt';
end
if 0 % rot(E)_x
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,rotEY_x)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '\nabla\times(E)_x';
end
if 0 % rot(E)_z
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,rotEY_z)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '\nabla\times(E)_z';
end
if 0
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,BX-rotEY_x)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'B_x(t1) + \Delta B_x';
end
if 0
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,0-rotEY_z)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'B_z(t1) + \Delta B_x';
end

if 0 % composite
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,EY1)
%hcb = colorbar('peer',hca);
%hcb.YLabel.String = 'B_z(t1)';
%hcb.YLabel.String = 'E_y';
colormap(hca,cmap)
if 1 % B
  hold(hca,'on')
  quiver(hca,Xq,Zq,BXq,BZq,'k','linewidth',1)
  hold(hca,'off')
end
if 1 % rotE
  hold(hca,'on')
  quiver(hca,Xq,Zq,-rotEYq_x,-rotEYq_z,'r','linewidth',1)
  hold(hca,'off')
end
legend(hca,{'E_y','B','-\nabla\times E'})
end
if 0 % composite
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,EY1)
%hcb = colorbar('peer',hca);
%hcb.YLabel.String = 'B_z(t1) + \Delta B_x';
%hcb.YLabel.String = 'E_y';
colormap(hca,cmap)
if 1 % B - rotE
  hold(hca,'on')
  quiver(hca,Xq,Zq,BXq-rotEYq_x,BZq-rotEYq_z,'k','linewidth',1)
  hold(hca,'off')
end
if 0 % rotE
  hold(hca,'on')
  quiver(hca,Xq,Zq,rotEYq_x,rotEYq_z,'r')
  hold(hca,'off')
end
legend(hca,{'E_y','B - {\Delta}t(\nabla\times E)'})
end

if 1 % composite, everything together
  hca = h(isub); isub = isub + 1;
  hE = pcolor(hca,X,Z,EY1);
  %hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'B_z(t1)';
  %hcb.YLabel.String = 'E_y';
  colormap(hca,cmap)
  if 1 % B-dt*rotX
    hold(hca,'on')
    hB2 = quiver(hca,Xq,Zq,BXq-rotEYq_x,BZq-rotEYq_z,0,'color',colors(3,:),'linewidth',1);
    hold(hca,'off')
  end
  if 1 % B
    hold(hca,'on')
    hB1 = quiver(hca,Xq,Zq,BXq,BZq,0,'color',0*colors(3,:),'linewidth',1);
    hold(hca,'off')
  end
  if 1 % rotE
    hold(hca,'on')
    hrotE = quiver(hca,Xq,Zq,-rotEYq_x,-rotEYq_z,0,'r','linewidth',1);
    hold(hca,'off')
  end
  legend([hE,hB1,hrotE,hB2],{'E_y','B','-\nabla\times E','B - {\Delta}t(\nabla\times E)'},'fontsize',16)
end

c_eval('shading(h(?),''flat'');',1:numel(h))
c_eval('axis(h(?),''equal'');',1:numel(h))
c_eval('h(?).Visible = ''off'';',1:numel(h))
hlinks = linkprop(h,{'XLim','YLim','CLim'});
hca.XLim = [-3.5 3.5];
hca.YLim = [-1.90 1.9];

%% Based on vector potential
B0 = 1;
L = 0.25;
E0 = 0.2;
Lx = 0.4;
Lz = 0.4;

t1 = 1;
t2 = 2.1;

Bx = @(z) B0*tanh(z/L);
Bx = @(z) B0*z./z;
Ay = @(z) B0*L*log(cosh(z/L));
%Ay = @(z) B0*L*z/L;
Ey = @(x,z,t) E0*exp(-x.^2/(Lx*t)^2-z.^2/(Lz*t)^2);



x = linspace(-3,3,500); dx = x(2) - x(1);
z = linspace(-1.6,1.6,400); dz = z(2) - z(1);
[X,Z] = ndgrid(x,z);


EY = Ey(X,Z,t1);
EY1 = Ey(X,Z,t1);
EY2 = Ey(X,Z,t2);
dA = EY;


% For quivers.
xq = linspace(-4,4,42*1); dxq = xq(2) - xq(1);
zq = linspace(-2,2,25*1); dzq = zq(2) - zq(1);
[Xq,Zq] = ndgrid(xq,zq);


BXq = Bx(Zq);
EYq = Ey(Xq,Zq,t1);
BZq = Zq*0;
%EY = Ey(X,Z,t2);

rotEYq_x = -[(EYq(:,2)-EYq(:,1))/dzq, (EYq(:,3:end)-EYq(:,1:end-2))/(2*dzq), (EYq(:,end)-EYq(:,end-1))/dzq];
rotEYq_z = [(EYq(2,:)-EYq(1,:))/dxq; (EYq(3:end,:)-EYq(1:end-2,:))/(2*dxq); (EYq(end,:)-EYq(end-1,:))/dxq];


% For quivers, 2.
xq2 = linspace(-4,4,8); dxq2 = xq2(2) - xq2(1);
zq2 = linspace(-2,2,7); dzq2 = zq2(2) - zq2(1);
[Xq2,Zq2] = ndgrid(xq2,zq2);

Xq2 = [0.4 -0.4];
Zq2 = [-1.17 1.17];

BXq2 = Bx(Zq2);
EYq2 = Ey(Xq2,Zq2,t1);
BZq2 = Zq2*0;

%% Plot

cmap = pic_colors('blue_red'); cmap = flipdim(cmap(1:floor(size(cmap,1)/2),:),1);
cmap = flipdim(pic_colors('thermal'),1);
cmap = flipdim(colormap('gray'),1);
colors = pic_colors('matlab');

alev = -1.5:0.1:1.5;

hca = subplot(1,1,1);
legs = {};

if 1 % A0  
  %hold(hca,'on')
  [hl_,hl] = contour(hca,X,Z,Ay(Z),alev,'k-','linewidth',1.0);
  %hold(hca,'off')  
  legs{end+1} = 'B';  
end
if 1 % Ey
  hold(hca,'on')
  hE = surf(hca,X,Z,X*0-1,EY2);
  shading(hca,'flat')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
  view([0 0 1])
  hold(hca,'off')
  legs{end+1} = 'E_y';
end
if 0 % B0
  hold(hca,'on')
  hB2 = quiver(hca,Xq2,Zq2,BXq2,BZq2,Sq,'color',colors(1,:),'linewidth',5);
  hold(hca,'off')
  legs{end+1} = '\Delta B';
end
if 1 % -rotX, ADD along stream linesaq
  hold(hca,'on')
  Sq = 0.5;
  lwidth = 1.5;
  hB2 = quiver(hca,Xq,Zq,-rotEYq_x,-rotEYq_z,Sq,'color',colors(3,:).^2,'linewidth',lwidth);
  hold(hca,'off')
  legs{end+1} = '\Delta{B}=-\Delta{t}\nabla\times E';
end
if 1 % A1
  hold(hca,'on')
  hl = contour(hca,X,Z,Ay(Z)+dA,alev,'k--','linewidth',1.0);  
  hold(hca,'off')
  legs{end+1} = 'B+\Delta{B}';
end






if 0 % B-dt*rotX, ADD along stream linesaq
  hold(hca,'on')
  hB2 = quiver(hca,Xq,Zq,BXq-rotEYq_x,BZq-rotEYq_z,Sq,'color','k','linewidth',lwidth);
  hold(hca,'off')
end

legend(hca,legs,'fontsize',14,'location','eastoutside','box','off')
colormap(hca,cmap)
hca.CLim = [0 E0*1.3];
hca.XLim = [-1.5 1.5];
axis equal
hca.XLim = [-1.5 1.5];
hca.YLim = [-1.1 1.1];
hca.Box = 'off';
axis off
hca.Position = [0.1300    0.1100    0.6023    0.8150];
%legend(hca,{'E_y','B','B-\int \nabla\times E dt','B','-\nabla\times E'},'fontsize',16,'location','eastoutside','box','off')

