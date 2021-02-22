%% Illustrate tube with helical lines

r = 0.5; % radius of tube
h = 2; % height of tube
axlim = 1;

% cylinder
nt = 100;
[X,Y,Z] = cylinder(r*0.999,nt); 
Z = Z*h;

% helical lines
nlines = 3;
nturns = 1.2;
ntheta = 200;
theta = linspace(0,nturns*360,ntheta);
z0 = 0;
z1 = h;
theta0 = linspace(0,360,nlines+1); theta0 = theta0(1:nlines);
x0 = r*cosd(theta0);
y0 = r*sind(theta0);

x = zeros(nlines,ntheta);
y = zeros(nlines,ntheta);
z = zeros(nlines,ntheta);
for iline = 1:nlines
  x(iline,:) = r*cosd(theta0(iline)+theta);
  y(iline,:) = r*sind(theta0(iline)+theta);
  z(iline,:) = linspace(0,h,ntheta);
end


% plot
hca = subplot(1,1,1);

h_patch = surf(hca,X,Y*1,Z); 
shading(hca,'flat')
h_patch.FaceAlpha = 0.8;

hold(hca,'on')
h_lines = plot3(x',y',z','linewidth',2,'color','k');
hold(hca,'off')

view(hca,[1 1 0.8])
axis(hca,'equal')
hca.XLim = axlim*[-1 1];
hca.YLim = axlim*[-1 1];
hca.ZLim = [0 h];
hca.Box = 'on';
hca.BoxStyle = 'full';

hca.XLabel.String = 'x';
hca.YLabel.String = 'y';
hca.ZLabel.String = 'z';

%% Rotating fluxtube
r = 0.2; % radius of tube
h = 2; % height of tube
axlim = 1;

% cylinder
nt = 100;
[X,Y,Z] = cylinder(r*0.999,nt); 
Z = Z*h;

% helical lines
nlines = 10;
nturns = 0.1;
ntheta = 200;

ntheta = 10;
theta0 = linspace(0,360,nlines+1); theta0 = theta0(1:nlines);
theta0_vec = linspace(0,2*360*0.2+360/40,40);


colors = pic_colors('matlab');
color1 = colors(1,:);
color2 = colors(2,:);
colortube = [0.9 0.9 0.95];
colortube = [0.95 0.95 1];

for it = 1:ntheta
  
  theta = linspace(0,nturns*360,ntheta);
  z0_ = 0;
  z1_ = 0.5*h;
  z2_ = h;
  
  theta1 = theta0;
  theta2 = theta0 + theta0_vec(it);
  x0 = r*cosd(theta0);
  y0 = r*sind(theta0);

  x1 = zeros(nlines,ntheta);
  y1 = zeros(nlines,ntheta);
  z1 = zeros(nlines,ntheta);
  for iline = 1:nlines
    x1(iline,:) = r*cosd(theta1(iline)+theta);
    y1(iline,:) = r*sind(theta1(iline)+theta);
    z1(iline,:) = linspace(z0_,z1_,ntheta);
  end
  x2 = zeros(nlines,ntheta);
  y2 = zeros(nlines,ntheta);
  z2 = zeros(nlines,ntheta);
  for iline = 1:nlines
    x2(iline,:) = r*cosd(theta2(iline)+theta);
    y2(iline,:) = r*sind(theta2(iline)+theta);
    z2(iline,:) = linspace(z1_,z2_,ntheta);
  end

  % plot
  hca = subplot(1,1,1);
  %hca.Visible = 'off';
  hca.Box = 'off';

  h_patch = surf(hca,X*0.98,Y*0.98,Z); 
  colormap(hca,colortube)
  shading(hca,'flat')
  h_patch.FaceAlpha = 0.6;

  hold(hca,'on')
  h_lines = plot3(hca,x1',y1',z1','linewidth',2,'color',color1);
  hpl1 = plot3(hca,x1(:,end),y1(:,end),z1(:,end),'color',color1,'Marker','o','MarkerFaceColor',color1,'linestyle','none')
  h_lines = plot3(hca,x2',y2',z2','linewidth',2,'color',color2);
  hpl2 = plot3(hca,x2(:,1),y2(:,1),z2(:,1),'color',color2,'Marker','o','MarkerFaceColor',color2,'linestyle','none');
  hold(hca,'off')
  
  hca.Box = 'off';
  

  view(hca,[1 1 0.8])
  axis(hca,'equal')
  hca.XLim = axlim*[-1 1];
  hca.YLim = axlim*[-1 1];
  hca.ZLim = [0 h];
  hca.Box = 'off';
  %hca.BoxStyle = 'full';
  hca.Visible = 'off';
  

  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
  drawnow
  pause(0.1)
end