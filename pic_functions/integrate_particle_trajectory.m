% integrate_particle_trajectory
doPlot = 1;
pic = df04;
xlim = pic.xi([1 pic.nx]); % di, should be used as conditions for when to stop itnegration
zlim = pic.zi([1 pic.nz]); % di

% Initial particle velocities and time when to start trajectory
disp('Preparing particles initial conditions.')
particleset = 4;
switch particleset
  case 1
    r0 = [100 0 0; 300 0 0]; % di
    t0 = [50, 50]; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 1;
    for iP = 1:nP      
      v0(iP,1) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vx(iSpecies);
      v0(iP,2) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vy(iSpecies);
      v0(iP,3) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vz(iSpecies);
    end
    v0 = [1 1 1; 1 1 1]; % vA
    t0 = [50, 50]; % wci-1
    T = [240, 240]; % wci-1
    
    m = [25, 1];
    q = [1, -1];
  case 2 % cold ions in inflow
    x_center = mean(pic.xi);
    x0 = x_center + (0:10);
    z0 = 3:5;
    [X0,Z0] = ndgrid(x0,z0);
    x0 = X0(:);
    z0 = Z0(:);
    r0 = [x0,x0*0,z0];
    t0 = zeros(nP,1); % wci-1
    T = t0 + 240; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 3; % cold ions from the north
    v0 = zeros(nP,3);
    for iP = 1:nP      
      v0(iP,1) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vx(iSpecies);
      v0(iP,2) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vy(iSpecies);
      v0(iP,3) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vz(iSpecies);
    end
    
    m = 25 + zeros(nP,1); % all ions
    q = 1 + zeros(nP,1);
  case 3 % hot ions from inflow
    x_center = mean(pic.xi);
    x0 = x_center + (0:5:10);
    z0 = 0:2:6;
    [X0,Z0] = ndgrid(x0,z0);
    x0 = X0(:);
    z0 = Z0(:);
    r0 = [x0,x0*0,z0];
    t0 = 4 + zeros(nP,1); % wci-1
    T = t0 + 240; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 1; % cold ions from the north
    v0 = zeros(nP,3);
    if 0
      for iP = 1:nP      
        v0(iP,1) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vx(iSpecies);
        v0(iP,2) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vy(iSpecies);
        v0(iP,3) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vz(iSpecies);
      end
    else
      for iP = 1:nP      
        v0(iP,1) = 0.1*randn(1,1);
        v0(iP,2) = 0.1*randn(1,1);
        v0(iP,3) = 0.1*randn(1,1);
      end
    end
    
    m = 25 + zeros(nP,1); % all ions
    q = 1 + zeros(nP,1);
  case 4 % hot electrons from inflow
    x_center = mean(pic.xi);
    x0 = x_center + (5:5:20);
    z0 = 0:2:6;
    [X0,Z0] = ndgrid(x0,z0);
    x0 = X0(:);
    z0 = Z0(:);
    r0 = [x0,x0*0,z0];
    nP = size(r0,1); % number of particles
    t0 = 4 + zeros(nP,1); % wci-1
    T = t0 + 240; % wci-1
    iSpecies = 1; % cold ions from the north
    v0 = zeros(nP,3);
    if 0
      for iP = 1:nP      
        v0(iP,1) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vx(iSpecies);
        v0(iP,2) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vy(iSpecies);
        v0(iP,3) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vz(iSpecies);
      end
    else
      for iP = 1:nP      
        v0(iP,1) = 0.2*randn(1,1);
        v0(iP,2) = 0.2*randn(1,1);
        v0(iP,3) = 0.2*randn(1,1);
      end
    end
    
    m = 1 + zeros(nP,1); % all ions
    q = -1 + zeros(nP,1);
end
%%
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

disp('Integrating trajectories.')
for iP = 1:nP  % one particle: 27s on office desktop
  tic  
  x_init = [r0(iP,:)'; v0(iP,:)']; % di, vA
  disp(sprintf('iP/nP = %g/%g, t0 = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',iP,nP,t0(iP),x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))
  
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


