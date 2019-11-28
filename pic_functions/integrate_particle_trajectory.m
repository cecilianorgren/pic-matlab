%% Initialize particles
doPlot = 1;
pic = df04;
xlim = pic.xi([1 pic.nx]); % di, should be used as conditions for when to stop itnegration
zlim = pic.zi([1 pic.nz]); % di

% Initial particle velocities and time when to start trajectory
disp('Preparing particles initial conditions.')
particleset = 2;
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
    t0 = 4 + zeros(nP,1); % wci-1
    T = t0 + 240; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 3; % cold ions from the north
    v0 = zeros(nP,3);
    if 0
      for iP = 1:nP      
        v0(iP,1) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vx(iSpecies);
        v0(iP,2) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vy(iSpecies);
        v0(iP,3) = pic.xlim(r0(iP,1)).zlim(r0(iP,3)).twcilim(t0(iP)).vz(iSpecies);
      end
    else
      for iP = 1:nP      
        v0(iP,1) = 0.02*randn(1,1);
        v0(iP,2) = 0.02*randn(1,1);
        v0(iP,3) = 0.02*randn(1,1);
      end
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
  case 5
    r0 = [220 0 0]; % di
    t0 = [50]; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 1;
    
    v0 = [0.1 0.1 0.1]; % vA
    t0 = [4]; % wci-1
    T = [240]; % wci-1
    
    m = [25];
    q = [1];
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
ttot = tic;
for iP = 1:nP  % one particle: 27s on office desktop
  tic  
  x_init = [r0(iP,:)'; v0(iP,:)']; % di, vA
  disp(sprintf('iP/nP = %g/%g, t0 = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',iP,nP,t0(iP),x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))
  
  % Integrate trajectory
  stopfunction = @(t,x,z) eom.box2d(t,x,z,xlim,zlim); % the stopfunction seems to require a lot of time, or not
  options = odeset('Events',stopfunction,'RelTol',1e-10);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
  options = odeset('InitialStep',0.05);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
  options = odeset('MaxStep',0.1);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
  options = odeset('RelTol',1e-14);
  
  EoM = @(ttt,xxx) eom_pic(ttt,xxx,pic,m(iP),q(iP));
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
  [t,x_sol] = ode45(EoM,[t0(iP) T(iP)],x_init,options);%,options); % 
  x_sol(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  
  x_sol_all{iP} = x_sol;
  
  toc
  if doPlot
    plot3(hca,x_sol(:,1),x_sol(:,2),x_sol(:,3))
    drawnow
  end
  
    if 0 % try smaller tolerance
      %%
    stopfunction = @(t,x,z) eom.box2d(t,x,z,xlim,zlim); % the stopfunction seems to require a lot of time, or not
    options = odeset('RelTol',1e-14);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

    EoM = @(ttt,xxx) eom_pic(ttt,xxx,pic,m(iP),q(iP));
    %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
    [t_,x_sol_] = ode45(EoM,[t0(iP) T(iP)],x_init,options); % ,options
    x_sol_(:,7) = t_; % x_sol = (x,y,z,vx,vy,vz,t)
    
    if doPlot
      plot3(hca,x_sol_(:,1),x_sol_(:,2),x_sol_(:,3))
      drawnow
    end
    end
  
end
if doPlot
  hold(hca,'off')
  hca.XLim = xlim;
  hca.ZLim = zlim;
end
xmax = 1;
toc(ttot)

%% Plot particles properties and forces as function of time

iP = 1;
nPanels = 6;
h = setup_subplots(nPanels,1);
isub = 1;
if 1 % xyz(t)  
  hca = h(isub); isub = isub + 1;
  plot(hca,x_sol(:,7),[x_sol(:,1)-x_sol(1,1) x_sol(:,2) x_sol(:,3)])
  hca.XLabel.String = 'twci';
  hca.YLabel.String = 'Position (d_i)';
  legend(hca,{'x-x(1)','y','z'})
end
if 1 % xyz(t)  
  hca = h(isub); isub = isub + 1;
  plot(hca,x_sol(:,7),[x_sol(:,3)])
  hca.XLabel.String = 'twci';
  hca.YLabel.String = 'Position (d_i)';
  legend(hca,{'z'})
end
if 1 % vxyz(t)  
  hca = h(isub); isub = isub + 1;
  plot(hca,x_sol(:,7),[x_sol(:,4) x_sol(:,5) x_sol(:,6)])
  hca.XLabel.String = 'twci';
  hca.YLabel.String = 'Velocity (v_A)';
  legend(hca,{'v_x','v_y','v_z'})
end
if 1 % Bxyz(t)  
  hca = h(isub); isub = isub + 1;
  plot(hca,x_sol(:,7),[Bx,By,Bz])
  hca.XLabel.String = 'twci';
  hca.YLabel.String = 'B (B_0)';
  legend(hca,{'B_x','B_y','B_z'})
end
if 1 % Exyz(t)  
  hca = h(isub); isub = isub + 1;
  plot(hca,x_sol(:,7),[Ex,Ey,Ez])
  hca.XLabel.String = 'twci';
  hca.YLabel.String = 'E (B_0v_A)';
  legend(hca,{'E_x','E_y','E_z'})
end
if 1 % ExB(t)  
  hca = h(isub); isub = isub + 1;
  plot(hca,x_sol(:,7),[Ey.*Bz-Ez.*By,Ez.*Bx-Ex.*Bz,Ex.*By-Ey.*Bz])
  hca.XLabel.String = 'twci';
  hca.YLabel.String = 'ExB (B_0^2v_A)';
  legend(hca,{'ExB_x','ExB_y','ExB_z'})
end

compact_panels(0.01)
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:nPanels)
%plot(x_sol(:,7),Ex,x_sol(:,7),Ey,x_sol(:,7),Ez,x_sol(:,7),Bx,x_sol(:,7),By,x_sol(:,7),Bz)

hlinks = linkprop(h,{'XLim'});

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


%% Compare forwards and backward integration for different error tolerances
% Initialize particles
doPlot = 1;
pic = df04;
xlim = pic.xi([1 pic.nx]); % di, should be used as conditions for when to stop itnegration
zlim = pic.zi([1 pic.nz]); % di

% Initial particle velocities and time when to start trajectory
disp('Preparing particles initial conditions.')
particleset = 5;
switch particleset
  case 5
    r0 = [220 0 0]; % di
    t0 = [50]; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 1;
    
    v0 = [0.1 0.1 0.1]; % vA
    t0 = [50]; % wci-1
    T = [240]; % wci-1
    
    m = [25];
    q = [1];
end

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
ttot = tic;
for iP = 1:nP  % one particle: 27s on office desktop
  options = odeset('RelTol',1e-18); % same for forward and backward
  %% Forward
  tic  
  x_init = [r0(iP,:)'; v0(iP,:)']; % di, vA
  disp(sprintf('iP/nP = %g/%g, t0 = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',iP,nP,t0(iP),x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))
  
  % Integrate trajectory  
  EoM = @(ttt,xxx) eom_pic(ttt,xxx,pic,m(iP),q(iP));  
  [t,x_sol_forw] = ode45(EoM,[t0(iP) T(iP)],x_init,options);
  x_sol_forw(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  toc
  
  %% Backward
  tic
  last_ind = find(not(isnan(x_sol_forw(:,1))),1,'last');
  x_init = tocolumn(x_sol_forw(last_ind,1:6)); % di, vA
  
  EoM = @(ttt,xxx) eom_pic_back(ttt,xxx,pic,m(iP),q(iP));
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
  [t,x_sol_back] = ode45(EoM,[x_sol_forw(last_ind,7) t0(iP)],x_init,options);%,options); % 
  x_sol_back(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  toc
  if doPlot
    %%
    plot3(hca,x_sol_forw(:,1),x_sol_forw(:,2),x_sol_forw(:,3),...
              x_sol_back(:,1),x_sol_back(:,2),x_sol_back(:,3),...
              x_sol_forw(1,1),x_sol_forw(1,2),x_sol_forw(1,3),'g*',...
              x_sol_forw(last_ind,1),x_sol_forw(last_ind,2),x_sol_forw(last_ind,3),'r*',...
              x_sol_back(1,1),x_sol_back(1,2),x_sol_back(1,3),'g+')
    drawnow
  end
  
end
if doPlot
  hold(hca,'off')
  hca.XLim = xlim;
  hca.ZLim = zlim;
end
xmax = 1;
toc(ttot)


%% I run doesnt go to T, restart it from last non nan instance
% Initialize particles
doPlot = 1;
pic = df04;
xlim = pic.xi([1 pic.nx]); % di, should be used as conditions for when to stop itnegration
zlim = pic.zi([1 pic.nz]); % di

% Initial particle velocities and time when to start trajectory
disp('Preparing particles initial conditions.')
particleset = 5;
switch particleset
  case 5
    r0 = [220 0 0]; % di
    t0 = [50]; % wci-1
    nP = size(r0,1); % number of particles
    iSpecies = 1;
    
    v0 = [0.1 0.1 0.1]; % vA
    t0 = [50]; % wci-1
    T = [240]; % wci-1
    
    m = [25];
    q = [1];
end

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
ttot = tic;
for iP = 1:nP  % one particle: 27s on office desktop
  options = odeset('RelTol',1e-14,'MaxStep',0.2); % same for forward and backward
  %% Forward
  tic  
  x_init = [r0(iP,:)'; v0(iP,:)']; % di, vA
  disp(sprintf('iP/nP = %g/%g, t0 = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',iP,nP,t0(iP),x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))
  
  % Integrate trajectory  
  EoM = @(ttt,xxx) eom_pic(ttt,xxx,pic,m(iP),q(iP));  
  [t,x_sol_1] = ode45(EoM,[t0(iP) T(iP)],x_init,options);
  x_sol_1(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  toc
  
  %% Continue
  tic
  last_ind = find(not(isnan(x_sol_1(:,1))),1,'last');
  x_init = tocolumn(x_sol_1(last_ind,1:6)); % di, vA
  
  EoM = @(ttt,xxx) eom_pic(ttt,xxx,pic,m(iP),q(iP));
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
  [t,x_sol_2] = ode45(EoM,[x_sol_1(last_ind,7) T(iP)],x_init,options);%,options); % 
  x_sol_2(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  toc
  %%
  if doPlot
    %%
    plot3(hca,x_sol_1(:,1),x_sol_1(:,2),x_sol_1(:,3),...
              x_sol_2(:,1),x_sol_2(:,2),x_sol_2(:,3),...
              x_sol_1(1,1),x_sol_1(1,2),x_sol_1(1,3),'g*',...
              x_sol_1(last_ind,1),x_sol_1(last_ind,2),x_sol_1(last_ind,3),'r*',...
              x_sol_2(1,1),x_sol_2(1,2),x_sol_2(1,3),'g+')
    drawnow
  end
  
end
if doPlot
  hold(hca,'off')
  hca.XLim = xlim;
  hca.ZLim = zlim;
end
xmax = 1;
toc(ttot)
