% demo for div_tensor and grad_scalar
localuser = datastore('local','user');
savedir_root = ['/Users/' localuser '/Research/PIC/df_cold_protons_1/'];

%% Load data
timestep = 8000;
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

%% Smoothing
structure_strs = {'pi1','pe1','pi2','pe2','nmvvi1','nmvve1','nmvvi2','nmvve2'};

%structure_strs = {'nmvvi1_smooth','nmvve1_smooth','nmvvi2_smooth','nmvve2_smooth'};
nstructures = numel(structure_strs);
for istructure = 1:nstructures  
  structure_str = structure_strs{istructure};
  structure_fields = {'xx','xy','xz','yy','yz','zz'};
  %structure_fields = {'xx','xy','xz','yy','yz','zz','scalar'};
  for ifield = 1:numel(structure_fields)
    eval(sprintf('%s_smooth.%s = smooth2(%s.%s,1);',structure_str,structure_fields{ifield},structure_str,structure_fields{ifield}))
  end
end

%% Make combines tensor T
T.xx = nmvvi1.xx + nmvve1.xx + nmvvi2.xx + nmvve2.xx + pi1.xx + pe1.xx + pi2.xx + pe2.xx - BB.xx + pB;
T.xy = nmvvi1.xy + nmvve1.xy + nmvvi2.xy + nmvve2.xy + pi1.xy + pe1.xy + pi2.xy + pe2.xy - BB.xy;
T.xz = nmvvi1.xz + nmvve1.xz + nmvvi2.xz + nmvve2.xz + pi1.xz + pe1.xz + pi2.xz + pe2.xz - BB.xz;
T.yy = nmvvi1.yy + nmvve1.yy + nmvvi2.yy + nmvve2.yy + pi1.yy + pe1.yy + pi2.yy + pe2.yy - BB.yy + pB;
T.yz = nmvvi1.yz + nmvve1.yz + nmvvi2.yz + nmvve2.yz + pi1.yz + pe1.yz + pi2.yz + pe2.yz - BB.yz;
T.zz = nmvvi1.zz + nmvve1.zz + nmvvi2.zz + nmvve2.zz + pi1.zz + pe1.zz + pi2.zz + pe2.zz - BB.zz + pB;
divT = div_tensor(x,z,T);

T_smooth.xx = nmvvi1_smooth.xx + nmvve1_smooth.xx + nmvvi2_smooth.xx + nmvve2_smooth.xx + pi1_smooth.xx + pe1_smooth.xx + pi2_smooth.xx + pe2_smooth.xx - BB.xx + pB;
T_smooth.xy = nmvvi1_smooth.xy + nmvve1_smooth.xy + nmvvi2_smooth.xy + nmvve2_smooth.xy + pi1_smooth.xy + pe1_smooth.xy + pi2_smooth.xy + pe2_smooth.xy - BB.xy;
T_smooth.xz = nmvvi1_smooth.xz + nmvve1_smooth.xz + nmvvi2_smooth.xz + nmvve2_smooth.xz + pi1_smooth.xz + pe1_smooth.xz + pi2_smooth.xz + pe2_smooth.xz - BB.xz;
T_smooth.yy = nmvvi1_smooth.yy + nmvve1_smooth.yy + nmvvi2_smooth.yy + nmvve2_smooth.yy + pi1_smooth.yy + pe1_smooth.yy + pi2_smooth.yy + pe2_smooth.yy - BB.yy + pB;
T_smooth.yz = nmvvi1_smooth.yz + nmvve1_smooth.yz + nmvvi2_smooth.yz + nmvve2_smooth.yz + pi1_smooth.yz + pe1_smooth.yz + pi2_smooth.yz + pe2_smooth.yz - BB.yz;
T_smooth.zz = nmvvi1_smooth.zz + nmvve1_smooth.zz + nmvvi2_smooth.zz + nmvve2_smooth.zz + pi1_smooth.zz + pe1_smooth.zz + pi2_smooth.zz + pe2_smooth.zz - BB.zz + pB;
divT_smooth = div_tensor(x,z,T_smooth);

%% Calculate gradients
% Ion pressure
gradpi1 = grad_scalar(x,z,pi1.scalar);
gradpi1_smooth = grad_scalar(x,z,pi1_smooth.scalar);
gradpi1_smooth_ = grad_scalar(x,z,smooth2(pi1.scalar,1));

%% Calculate divergences
% Pressure
c_eval('divpe? = div_tensor(x,z,pe?);',1:2)
c_eval('divpe?_smooth = div_tensor(x,z,pe?_smooth);',1:2)
c_eval('divpi? = div_tensor(x,z,pi?);',1:2)
c_eval('divpi?_smooth = div_tensor(x,z,pi?_smooth);',1:2)

% Dynamic pressure
c_eval('divnmvvi? = div_tensor(x,z,nmvvi?);',1:2)
c_eval('divnmvvi?_smooth = div_tensor(x,z,nmvvi?_smooth);',1:2)
c_eval('divnmvve? = div_tensor(x,z,nmvve?);',1:2)
c_eval('divnmvve?_smooth = div_tensor(x,z,nmvve?_smooth);',1:2)

%c_eval('divnmvvi?_smooth_smooth = div_tensor(x,z,nmvvi?_smooth_smooth);',1:2)
%c_eval('divnmvve?_smooth_smooth = div_tensor(x,z,nmvve?_smooth_smooth);',1:2)

% Magnetic field
divBB = div_tensor(x,z,BB);
gradpB = grad_scalar(x,z,pB);

%% Plot divergence of T components
xlim = [0 60];%[x(1) x(end)] + 150*[1 -1];
zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

% Flux function
doA = 1;
cA = 0*[0.8 0.8 0.8];
nA = 20;
nA = [0:-1:min(A(:))];
ipxA = ipx1:20:ipx2;
ipzA = ipz1:20:ipz2;


% Variables to plot
varstrs = {'divpi1_smooth.x','divpi1_smooth.y','divpi1_smooth.z',...
  'divpe1_smooth.x','divpe1_smooth.y','divpe1_smooth.z',...
  'divpi2_smooth.x','divpi2_smooth.y','divpi2_smooth.z',...
  'divpe2_smooth.x','divpe2_smooth.y','divpe2_smooth.z',...
  ...%'divnmvvi1_smooth.x','divnmvvi1_smooth.y','divnmvvi1_smooth.z',...
  ...%'divnmvvi2_smooth.x','divnmvvi2_smooth.y','divnmvvi2_smooth.z',...
  ...%'gradpB.x','gradpB.y','gradpB.z',...
  'divBB.x','divBB.y','divBB.z'...
  };
varstrs = {'divpi1_smooth.x','divpi1_smooth.y','divpi1_smooth.z',...
  'divpe1_smooth.x','divpe1_smooth.y','divpe1_smooth.z',...
  'divpi2_smooth.x','divpi2_smooth.y','divpi2_smooth.z',...
  'divpe2_smooth.x','divpe2_smooth.y','divpe2_smooth.z',...
  'divpi1_smooth.x+divpe1_smooth.x+divpi2_smooth.x+divpe2_smooth.x',...
  'divpi1_smooth.y+divpe1_smooth.y+divpi2_smooth.y+divpe2_smooth.y',...
  'divpi1_smooth.z+divpe1_smooth.z+divpi2_smooth.z+divpe2_smooth.z'...
  };
varstrs = {'divnmvvi1_smooth.x','divnmvvi1_smooth.y','divnmvvi1_smooth.z',...
  'divnmvve1_smooth.x','divnmvve1_smooth.y','divnmvve1_smooth.z',...
  'divnmvvi2_smooth.x','divnmvvi2_smooth.y','divnmvvi2_smooth.z',...
  'divnmvve2_smooth.x','divnmvve2_smooth.y','divnmvve2_smooth.z',...
  'divnmvvi1_smooth.x+divnmvve1_smooth.x+divnmvvi2_smooth.x+divnmvve2_smooth.x',...
  'divnmvvi1_smooth.y+divnmvve1_smooth.y+divnmvvi2_smooth.y+divnmvve2_smooth.y',...
  'divnmvvi1_smooth.z+divnmvve1_smooth.z+divnmvvi2_smooth.z+divnmvve2_smooth.z'...
  };
varstrs = {'nmvvi1_smooth.xx','nmvvi1_smooth.xy','nmvvi1_smooth.xz',...
  'nmvvi1_smooth.yy','nmvvi1_smooth.yz','nmvvi1_smooth.zz',...
  'divnmvvi1_smooth.x','divnmvvi1_smooth.y','divnmvvi1_smooth.z',...
  };
varstrs = {'nmvvi1.xx','nmvvi1.xy','nmvvi1.xz',...
  'nmvvi1.yy','nmvvi1.yz','nmvvi1.zz',...
  'divnmvvi1.x','divnmvvi1.y','divnmvvi1.z',...
  };
varstrs = {'nmvvi1.xx','nmvvi1.xy','nmvvi1.xz',...
  'nmvvi1.yy','nmvvi1.yz','nmvvi1.zz',...
  'nmvvi1_smooth.xx','nmvvi1_smooth.xy','nmvvi1_smooth.xz',...
  'nmvvi1_smooth.yy','nmvvi1_smooth.yz','nmvvi1_smooth.zz'
  };
varstrs = {'divnmvvi1.x','divnmvvi1.y','divnmvvi1.z',...
  'divnmvvi1_smooth.x','divnmvvi1_smooth.y','divnmvvi1_smooth.z',...
  'divnmvvi1.x-divnmvvi1_smooth.x','divnmvvi1.y-divnmvvi1_smooth.y','divnmvvi1.z-divnmvvi1_smooth.z'
  };
varstrs = {'divnmvvi1.x','divnmvvi1.y','divnmvvi1.z',...
  'divnmvvi1_smooth.x','divnmvvi1_smooth.y','divnmvvi1_smooth.z',...
  'divnmvvi1_smooth_smooth.x','divnmvvi1_smooth_smooth.y','divnmvvi1_smooth_smooth.z'%,...
  ...%'divnmvvi1.x-divnmvvi1_smooth.x','divnmvvi1.y-divnmvvi1_smooth.y','divnmvvi1.z-divnmvvi1_smooth.z'
  };
varstrs = {...
  'divpi1_smooth.x','divpi1_smooth.y','divpi1_smooth.z',...
  'divnmvvi1_smooth_smooth.x','divnmvvi1_smooth_smooth.y','divnmvvi1_smooth_smooth.z',...  
  'gradpB.x','gradpB.y','gradpB.z',...
  '-divBB.x','-divBB.y','-divBB.z'...
  };
varstrs = {...
  'divpi1_smooth.x','divpi1_smooth.y','divpi1_smooth.z',...
  'divnmvvi1_smooth_smooth.x','divnmvvi1_smooth_smooth.y','divnmvvi1_smooth_smooth.z',...  
  'divpi2_smooth.x','divpi2_smooth.y','divpi2_smooth.z',...
  'divnmvvi2_smooth_smooth.x','divnmvvi2_smooth_smooth.y','divnmvvi2_smooth_smooth.z',...  
  'gradpB.x','gradpB.y','gradpB.z',...
  '-divBB.x','-divBB.y','-divBB.z'...
  };
varstrs = {...  
  %'B.x','B.y','B.z',...
  'vi1.par','E.par','E_smooth.par',...
  'vExB.x','vExB.y','vExB.z',...
  '-divBB.x','-divBB.y','-divBB.z'...
  };
varstrs = {...      
  ...%'divT.x','divT.y','divT.z',...
  'T_smooth.xx','T_smooth.xy','T_smooth.xz'...
  'T_smooth.yy','T_smooth.yz','T_smooth.zz'...
  };
varstrs = {...      
  ...%'divT.x','divT.y','divT.z',...
  'divT_smooth.x','divT_smooth.y','divT_smooth.z'...
  };
varstrs = {'divpi1_smooth.x','divpi1_smooth.y','divpi1_smooth.z'};
%varstrs = {'divpe1_smooth.x','divpe1_smooth.y','divpe1_smooth.z'};
%varstrs = {'divpi2_smooth.x','divpi2_smooth.y','divpi2_smooth.z'};
%varstrs = {'divpe2_smooth.x','divpe2_smooth.y','divpe2_smooth.z'};

varstrs = {'divnmvvi1_smooth.x','divnmvvi1_smooth.y','divnmvvi1_smooth.z'};
varstrs = {'divnmvve1_smooth.x','divnmvve1_smooth.y','divnmvve1_smooth.z'};
varstrs = {'divnmvvi2_smooth.x','divnmvvi2_smooth.y','divnmvvi2_smooth.z'};
varstrs = {'divnmvve2_smooth.x','divnmvve2_smooth.y','divnmvve2_smooth.z'};
varstrs = {'gradpB.x','gradpB.y','gradpB.z'};
varstrs = {'-divBB.x','-divBB.y','-divBB.z'};
varstrs = {'divT_smooth.x','divT_smooth.y','divT_smooth.z'};
%varstrs = {'divT.x','divT.y','divT.z'};

% Quivers
doQ = 1;
nQx = 80;
nQz = 20;
[Z,X] = meshgrid(z,x);
ipxQ = fix(linspace(ipx1,ipx2,nQx));
ipzQ = fix(linspace(ipz1,ipz2,nQz));
[dataQx,dataQz] = meshgrid(ipxQ,ipzQ);
ipXQ = dataQx; ipZQ = dataQz;
dataQ = divT_smooth;
%dataQ.x = -divBB.x;
%dataQ.y = -divBB.y;
%dataQ.z = -divBB.z;
maxQ = 1;
dataQ.abs = sqrt(dataQ.x.^2 + dataQ.z.^2);
dataQ.x(dataQ.abs>maxQ) = NaN;
dataQ.y(dataQ.abs>maxQ) = NaN;
dataQ.z(dataQ.abs>maxQ) = NaN;


% Set up figure
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

nvars = numel(varstrs);
for ivar = 1:nvars
  hca = h(isub); isub = isub + 1;
  varstr = varstrs{ivar};
  variable = eval(varstr);
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.XLabel.String = 'x (di)';
  hca.YLabel.String = 'z (di)';
  hcb = colorbar('peer',hca);
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';
  colormap(hca,pic_colors('blue_red'))
  %colormap(hca,pic_colors('candy2'))
  hca.XDir = 'reverse';
  hca.YDir = 'normal';
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    hold(hca,'off')  
  end
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hquiv.LineWidth = 1;
    hquiv.Color = 0*[0.9290    0.6940    0.1250];
    hold(hca,'off')  
  end
end
hlink = linkprop(h,{'CLim','XLim','YLim'});
c_eval('axis(h(?),''equal'')',1:npanels)
%hlink = linkprop(h,{'XLim','YLim'});
h(1).CLim = [-0.3 0.3];

% %% Step 1, compare div_tensor(tens_P) and grad_scalar(scal_P)
%% Plot pi1.scalar and pi1_smooth.scalar
xlim = [0 60];%[x(1) x(end)] + 150*[1 -1];
zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;


varstrs = {'pi1.scalar','pi1_smooth.scalar','pi1.scalar-pi1_smooth.scalar'};
nvars = numel(varstrs);
for ivar = 1:nvars
  hca = h(isub); isub = isub + 1;
  varstr = varstrs{ivar};
  variable = eval(varstr);
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.XLabel.String = 'x (di)';
  hca.YLabel.String = 'z (di)';
  hcb = colorbar('peer',hca);
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';
  colormap(hca,pic_colors('blue_red'))
end
hlink = linkprop(h(1:2),{'CLim','XLim','YLim'});
h(1).CLim = [-1 1];

%% Plot divergence and gradient
xlim = [0 60];%[x(1) x(end)] + 150*[1 -1];
zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;


varstrs = {'gradpi1.x','gradpi1.y','gradpi1.z',...
  'gradpi1_smooth_.x','gradpi1_smooth_.y','gradpi1_smooth_.z',...
  'divpi1.x','divpi1.y','divpi1.z',...
  'divpi1_smooth.x','divpi1_smooth.y','divpi1_smooth.z'...
  ...%'pi1.scalar','pi1_smooth.scalar'
  };
varstrs = {'divBB.x','divBB.y','divBB.z'};
nvars = numel(varstrs);
for ivar = 1:nvars
  hca = h(isub); isub = isub + 1;
  varstr = varstrs{ivar};
  variable = eval(varstr);
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.XLabel.String = 'x (di)';
  hca.YLabel.String = 'z (di)';
  hcb = colorbar('peer',hca);
  hca.Title.String = varstr;
  hca.Title.Interpreter = 'none';
  colormap(hca,pic_colors('blue_red'))
  hca.XDir = 'reverse';
  hca.YDir = 'normal';
end
hlink = linkprop(h,{'CLim','XLim','YLim'});
h(1).CLim = [-1 1];


%%
T.xx = nmvvi1.xx + nmvve1.xx + nmvvi2.xx + nmvve2.xx + pi1.xx + pe1.xx + pi2.xx + pe2.xx + BB.xx + pB;
T.xy = nmvvi1.xy + nmvve1.xy + nmvvi2.xy + nmvve2.xy + pi1.xy + pe1.xy + pi2.xy + pe2.xy + BB.xy;
T.xz = nmvvi1.xz + nmvve1.xz + nmvvi2.xz + nmvve2.xz + pi1.xz + pe1.xz + pi2.xz + pe2.xz + BB.xz;
T.yy = nmvvi1.yy + nmvve1.yy + nmvvi2.yy + nmvve2.yy + pi1.yy + pe1.yy + pi2.yy + pe2.yy + BB.yy + pB;
T.yz = nmvvi1.yz + nmvve1.yz + nmvvi2.yz + nmvve2.yz + pi1.yz + pe1.yz + pi2.yz + pe2.yz + BB.yz;
T.zz = nmvvi1.zz + nmvve1.zz + nmvvi2.zz + nmvve2.zz + pi1.zz + pe1.zz + pi2.zz + pe2.zz + BB.zz + pB;
divT = div_tensor(x,z,T);