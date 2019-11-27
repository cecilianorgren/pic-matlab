% integrate_particle_trajectory
doPlot = 1;
pic = df04;
xlim = pic.xi([1 pic.nx]); % di, should be used as conditions for when to stop itnegration
zlim = pic.zi([1 pic.nz]); % di

% Initial particle velocities and time when to start trajectory
r0 = [100 0 0; 300 0 0]; % di
v0 = [1 1 1; 1 1 1]; % vA
t0 = [50, 50]; % wci-1
T = [200, 200]; % wci-1
nP = size(r0,1); % number of particles
m = [25, 1];
q = [1, -1];

x_sol_all = cell(nP,1);

if doPlot
  hca = subplot(1,1,1);  
  hca.XLim = xlim;
  hca.ZLim = zlim;
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'y/d_i';
  hca.ZLabel.String = 'z/d_i';
  hold(hca,'on')
end
for iP = 1:nP  % one particle: 27s on office desktop
  tic
  x_init = [r0(iP,:)'; v0(iP,:)']; % di, vA
  
  % Integrate trajectory
  stopfunction = @(t,x,z) eom.box2d(t,x,z,xlim,zlim); % the stopfunction seems to require a lot of time, or not
  options = odeset('Events',stopfunction,'RelTol',1e-10);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

  EoM = @(ttt,xxx) eom_pic(ttt,xxx,pic,m(iP),q(iP));
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
  [t,x_sol] = ode45(EoM,[t0(iP) T(iP)],x_init); % ,options
  x_sol(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  
  %x = x_sol(:,1);
  %y = x_sol(:,2);
  %z = x_sol(:,3);
  %vx = x_sol(:,4);
  %vy = x_sol(:,5);
  %vz = x_sol(:,6);
  
  
  
  x_sol_all{iP} = x_sol;
  toc
  if doPlot
    plot3(hca,x_sol(:,1),x_sol(:,2),x_sol(:,3))
    drawnow
  end
end
if doPlot
  hold(hca,'off')
  hca.XLim = xlim;
  hca.ZLim = zlim;
end
xmax = 1;
%% Plot particles on top of fields

h = setup_subplots(1,1,1);
isub = 1;
for it = 1:pic.length
  hca = h(isub); isub = isub + 1;
  xmax = max(x_sol_all{iP}(:,1));
  xmin = min(x_sol_all{iP}(:,1));
  zmax = max(x_sol_all{iP}(:,3));
  zmin = min(x_sol_all{iP}(:,3));
  pic = df04.xlim([xmin xmax]+[-2 2]).zlim([zmin zmax]+[-2 2]);
  imagesc(hca,pic.xi,pic,pic.Ey')
  hold(hca,'on')
  for iP = 1:nP
    plot3(hca,x_sol_all{iP}(:,1),x_sol_all{iP}(:,2),x_sol_all{iP}(:,3))
  end
  hold(hca,'off')
  pause(0.1)
end


