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
