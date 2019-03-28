
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/';

%% Load data
timestep = 08000;
txtfile = sprintf('/Users/cecilia/Data/PIC/data/fields-%05.0f.dat',timestep); % michael's perturbation
%txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons_1/data/fields-%05.0f.dat',timestep); % michael's perturbation
txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/fields-%05.0f.dat',timestep); % michael's perturbation

%timestep = 055;
%txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons_2/data/fields-%05.0f.dat',timestep); % try-out with larger perturbation


%tic; [time,r,e,b,n1,ne,ve,vi,je,ji,pe,pi,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q] = read_fields_ieie(txtfile); toc
tic; [x,z,E,B,...s
  ni1,ne1,ni2,ne2,...
  vi1,ve1,vi2,ve2,...
  ji1,je1,ji2,je2,...
  pi1,pe1,pi2,pe2,...
  ti1,te1,ti2,te2,...
  dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile); toc
x0 = mean(x); x = x-x0;

% Calculate auxillary quantities
A = vector_potential(x,z,B.x,B.z); % vector potential
%phi = scalar_potential(x,z,E.x,E.z); % vector potential
[saddle_locations,saddle_values] = saddle(A);
%[A_volume,A_map] = fluxtube_volume(A,-30:1:1);
%[saddle_locations,saddle_values] = fluxtube_volume(A);

pic_calc_script

r1 = b; % magnetic field unit vector
r2 = cross_product(r1.x,r1.y,r1.z,0,1,0);
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2.x = r2.x./r2.abs;
r2.y = r2.y./r2.abs;
r2.z = r2.z./r2.abs;
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2 = cross_product(r2.x,r2.y,r2.z,r1.x,r1.y,r1.z);
r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);
r3.abs = sqrt(r3.x.^2 + r3.y.^2 + r3.z.^2);

tic;te1_fac = rotate_tens(te1,r1,r2,r3); toc
tic;te2_fac = rotate_tens(te2,r1,r2,r3); toc
tic;ti1_fac = rotate_tens(ti1,r1,r2,r3); toc
tic;ti2_fac = rotate_tens(ti2,r1,r2,r3); toc
te1_perp = 0.5*(te1.yy + te1.zz);
te1_par = te1.xx;
te2_perp = 0.5*(te2.yy + te2.zz);
te2_par = te2.xx;
ti1_perp = 0.5*(ti1.yy + ti1.zz);
ti1_par = ti1.xx;
ti2_perp = 0.5*(ti2.yy + ti2.zz);
ti2_par = ti2.xx;

% Stream functions
c_eval('Se?.xz = vector_potential(x,z,ve?.x,ve?.z);',1:2) % stream function
c_eval('Si?.xz = vector_potential(x,z,vi?.x,vi?.z);',1:2) % stream function

jtot.x = ji1.x + ji2.x - je1.x - je2.x;
jtot.y = ji1.y + ji2.y - je1.y - je2.y;
jtot.z = ji1.z + ji2.z - je1.z - je2.z;


%% Plot, define variable in cell array
% Define what variables to plot
varstrs = {'ve1.x','ve2.x','ve1.z','ve2.z','ve1.par','ve2.par','ExB.x','ExB.z','-ve1xB.x','-ve2xB.x','-ve1xB.z','-ve2xB.z','E.x','E.z'};
varstrs = {'ve1.x','ve2.x','B.z','E.z','-ve1xB.x','A'};
varstrs = {'ve1.perp.x','ve2.perp.x','vi1.perp.x','vi2.perp.x','ExB.x'};  
varstrs = {'ve1.perp.z','ve2.perp.z','vi1.perp.z','vi2.perp.z','ExB.z'};
varstrs = {'ne1','ne2','ni1','ni2','te2.scalar','ti2.scalar','pe2.scalar','pi2.scalar'};
varstrs = {'ne1','ne2','ni1','ni2'};

varstrs = {'Ute1','Ute2','Uti1','Uti2','Uke1','Uke2','Uki1','Uki2'}; clim = 12*[-1 1];
varstrs = {'pe1.xx','pe1.xy','pe1.yy','pe1.xz','pe1.zz','pe1.yz'}; clim = 0.25*[-1 1];
varstrs = {'ne1','ne2','ni1','ni2'}; clim = 2*[-1 1];
varstrs = {'pe2.xx','pe2.xy','pe2.yy','pe2.xz','pe2.zz','pe2.yz'}; clim = 0.25*[-1 1];
varstrs = {'pi2.xx','pi2.xy','pi2.yy','pi2.xz','pi2.zz','pi2.yz'}; clim = 0.3*[-1 1];
varstrs = {'ve1.x','ve1.y','ve1.z','ExB.x','ExB.y','ExB.z','-ve1xB.x','-ve1xB.y','-ve1xB.z','-E.x','-E.y','-E.z','-ve1xB.x-E.x','-ve1xB.y-E.y','-ve1xB.z-E.z'}; clim = [-1 1];
varstrs = {'vi1.x','vi1.y','vi1.z','ExB.x','ExB.y','ExB.z','vi1xB.x','vi1xB.y','vi1xB.z','E.x','E.y','E.z','vi1xB.x+E.x','vi1xB.y+E.y','vi1xB.z+E.z'}; clim = [-1 1];
varstrs = {'pi1.xx','pi1.xy','pi1.yy','pi1.xz','pi1.zz','pi1.yz','gradpi1.x','gradpi1.y','gradpi1.z'}; clim = 0.3*[-1 1];
varstrs = {'pi1.xx','pi1.zz','gradpi1.x','gradpi1.z','gradpi1_smooth.x','gradpi1_smooth.z','gradpi1.abs','gradpi1_smooth.abs'}; clim = 0.3*[-1 1];
varstrs = {'-ne1.*ve1xB.x','-ne1.*ve1xB.y','-ne1.*ve1xB.z','-ne1.*E.x','-ne1.*E.y','-ne1.*E.z','-ne1.*(ve1xB.x+E.x)','-ne1.*(ve1xB.y+E.y)','-ne1.*(ve1xB.z+E.z)','-gradpe1_smooth.x','-gradpe1_smooth.y','-gradpe1_smooth.z'}; clim = 0.5*[-1 1];
varstrs = {'-ne2.*ve2xB.x','-ne2.*ve2xB.y','-ne2.*ve2xB.z','-ne2.*E.x','-ne2.*E.y','-ne2.*E.z',...
           '-ne2.*(ve2xB.x+E.x)','-ne2.*(ve2xB.y+E.y)','-ne2.*(ve2xB.z+E.z)',...
           '-gradpe2_smooth.x','-gradpe2_smooth.y','-gradpe2_smooth.z',...
           '-ne2.*(ve2xB.x+E.x)-gradpe2_smooth.x','-ne2.*(ve2xB.y+E.y)-gradpe2_smooth.y','-ne2.*(ve2xB.z+E.z)-gradpe2_smooth.z'...
           }; clim = 0.5*[-1 1];
varstrs = {'ve1.x','ve2.x','vi1.x','vi2.x'}; clim = 3*[-1 1];
varstrs = {'te1.scalar','te2.scalar','te2.scalar./te1.scalar','ti1.scalar','ti2.scalar','ti2.scalar./ti1.scalar'}; clim = 0.8*[-1 1];
varstrs = {'te1.scalar','te2.scalar','te1.scalar./te2.scalar','ne1','ne2'}; clim = 3*[-1 1];
varstrs = {'ve1.par','ve1.perp.x','ve1.perp.y','ve1.perp.z','ve1.par./sqrt(ve1.perp.x.^2+ve1.perp.y.^2+ve1.perp.z.^2)'}; clim = 3*[-1 1];
%varstrs = {'vi1.par','vi1.perp.x','vi1.perp.y','vi1.perp.z','vi1.par./sqrt(vi1.perp.x.^2+vi1.perp.y.^2+vi1.perp.z.^2)'}; clim = 3*[-1 1];
varstrs = {'ne1','ne2'}; clim = 3*[-1 1];
varstrs = {'(-ve1xB.y+ve2xB.y)','(-vi1xB.y+vi2xB.y)','vi1.y','vi2.y'}; clim = 0.1*[-1 1];
varstrs = {'(-vi1xB.y+vi2xB.y)','vi1.y','vi2.y'}; clim = [];0.2*[-1 1];
varstrs = {'vi2xB.y','vi2xB.y_zx','vi2xB.y_xz'}; clim = [];0.2*[-1 1];
varstrs = {'vi1xB.y','vi2xB.y','vi1xB.y_zx','vi2xB.y_zx','vi1xB.y_xz','vi2xB.y_xz'}; clim = 0.5*[-1 1];
varstrs = {'vi1xB.y','vi2xB.y','-vi1xB.y+vi2xB.y','vi1xB.y_zx','vi2xB.y_zx','-vi1xB.y_zx+vi2xB.y_zx','vi1xB.y_xz','vi2xB.y_xz','-vi1xB.y_xz+vi2xB.y_xz'}; clim = 0.2*[-1 1];
%varstrs = {'-vi1xB.y+vi2xB.y','-vi1xB.y_zx+vi2xB.y_zx','-vi1xB.y_xz+vi2xB.y_xz','vi1.x','vi1.y','vi1.z','vi2.x','vi2.y','vi2.z'}; clim = 0.2*[-1 1];
%varstrs = {'vte1','vte2','vti1','vti2'}; clim = [];0.2*[-1 1];
%varstrs = {'wce1','wce2','wci1','wci2','vte1','vte2','vti1','vti2','re1','re2','ri1','ri2'}; clim = 10*[-1 1];0.2*[-1 1];
%varstrs = {'vte1','vte2','vti1','vti2','re1','re2','ri1','ri2'}; clim = 3*[-1 1];0.2*[-1 1];
varstrs = {'pi1.scalar','gradpi1.x','gradpi1.z','vExB.x','vExB.y','vExB.z'}; clim = [];
varstrs = {'vExB.x_yz','vExB.x_zy','vExB.y_zx','vExB.y_xz','vExB.z_xy','vExB.z_yx'}; clim = [];
varstrs = {'vExB.x','vExB.y','vExB.z','vExB.x_yz','vExB.x_zy','vExB.y_zx','vExB.y_xz','vExB.z_xy','vExB.z_yx'}; clim = [-2 2];
varstrs = {'vExB.x','vExB.x_yz','vExB.x_zy','vExB.y','vExB.y_zx','vExB.y_xz','vExB.z','vExB.z_xy','vExB.z_yx'}; clim = [-2 2];

varstrs = {'vExB.x','vExB.y','vExB.z','vDi1.x','vDi1.y','vDi1.z'}; clim = [];
varstrs = {'vExB.x','vExB.y','vExB.z','vDi1.x','vDi1.y','vDi1.z','vDi2.x','vDi2.y','vDi2.z','vDe1.x','vDe1.y','vDe1.z','vDe2.x','vDe2.y','vDe2.z'}; clim = [-1.1 1.1];
varstrs = {'vi1.y','vDi1.y','vExB.y','vDi1.y+vExB.y','vi2.y','vDi2.y','vExB.y','vDi2.y+vExB.y','ve1.y','vDe1.y','vExB.y','vDe1.y+vExB.y','ve2.y','vDe2.y','vExB.y','vDe2.y+vExB.y'}; clim = [-1.1 1.1];
varstrs = {'ve.x','ve.y','ve.z','vi.x','vi.y','vi.z'}; clim = 3*[-1 1];
varstrs = {'vi1.x','ve1.x','vi2.x','ve2.x','vExB.x','vi1.x-vDF','ve1.x-vDF','vi2.x-vDF','ve2.x-vDF','vExB.x-vDF','vi1.x-vExB.x','ve1.x-vExB.x','vi2.x-vExB.x','ve2.x-vExB.x','vExB.x-vExB.x'}; clim = 2*[-1 1];
varstrs = {'vi1.z','ve1.z','vi2.z','ve2.z','vExB.z','vi1.x-vDF','ve1.x-vDF','vi2.x-vDF','ve2.x-vDF','vExB.x-vDF','vi1.x-vExB.x','ve1.x-vExB.x','vi2.x-vExB.x','ve2.x-vExB.x','vExB.x-vExB.x'}; clim = 2*[-1 1];
varstrs = {'B.y'}; clim = 0.5*[-1 1];
varstrs = {'Uti1-Uki1','Uke1-Ute1','Uki2-Uti2','Uke2-Ute2+0.02'}; clim = [-0.3 0.3];
%varstrs = {'ji1.x','je1.x','ji1.x-je1.x','ji2.x','je2.x','ji2.x-je2.x','ji1.x-je1.x+ji2.x-je2.x'}; clim = 1.2*[-1 1];
%varstrs = {'ne1','ne2','ne','ni1','ni2','ni'}; clim = [];
varstrs = {'T_smooth.xx','T_smooth.xy','T_smooth.xz','T_smooth.yy','T_smooth.yz','T_smooth.zz'}; clim = 1.5*[-1 1];
varstrs = {'nmvvi1_smooth.xx ','nmvve1_smooth.xx','nmvvi2_smooth.xx','nmvve2_smooth.xx','pi1_smooth.xx','pe1_smooth.xx','pi2_smooth.xx','pe2_smooth.xx','-BB.xx','pB'}; clim = 1.5*[-1 1];
%varstrs = {'nmvvi1_smooth.xy ','nmvve1_smooth.xy','nmvvi2_smooth.xy','nmvve2_smooth.xy','pi1_smooth.xy','pe1_smooth.xy','pi2_smooth.xy','pe2_smooth.xy','-BB.xy'}; clim = 1.5*[-1 1];
%varstrs = {'nmvvi1_smooth.xz ','nmvve1_smooth.xz','nmvvi2_smooth.xz','nmvve2_smooth.xz','pi1_smooth.xz','pe1_smooth.xz','pi2_smooth.xz','pe2_smooth.xz','-BB.xz'}; clim = 1.5*[-1 1];
%varstrs = {'nmvvi1_smooth.yy ','nmvve1_smooth.yy','nmvvi2_smooth.yy','nmvve2_smooth.yy','pi1_smooth.yy','pe1_smooth.yy','pi2_smooth.yy','pe2_smooth.yy','-BB.yy','pB'}; clim = 1.5*[-1 1];
%varstrs = {'nmvvi1_smooth.yz ','nmvve1_smooth.yz','nmvvi2_smooth.yz','nmvve2_smooth.yz','pi1_smooth.yz','pe1_smooth.yz','pi2_smooth.yz','pe2_smooth.yz','-BB.yz'}; clim = 1.5*[-1 1];
%varstrs = {'nmvvi1_smooth.zz ','nmvve1_smooth.zz','nmvvi2_smooth.zz','nmvve2_smooth.zz','pi1_smooth.zz','pe1_smooth.zz','pi2_smooth.zz','pe2_smooth.zz','-BB.zz','pB'}; clim = 1.5*[-1 1];

nvars = numel(varstrs);

%xlim = torow(x([1 end])) + [100 -100];
%zlim = [-15 15];z([1 end]);

% Initialize figure
npanels = nvars;
nrows = 4;
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end
%h = setup_subplots(nrows,ncols,'horizontal');
h = setup_subplots(nrows,ncols,'vertical');
linkaxes(h);

xlim = [x(1) x(end)]; xlim = [0 70];%[x(1) x(end)] + 150*[1 -1];
zlim = [z(1) z(end)]; zlim = 10*[-1 1];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:2:ipx2;
ipz = ipz1:2:ipz2;
    
% Flux function
doA = 1;
cA = 0*[0.8 0.8 0.8];
nA = 20;
nA = [0:-1:min(A(:))];
ipxA = ipx1:20:ipx2;
ipzA = ipz1:20:ipz2;

%sepA = A(find(B.abs(:)==min(B.abs(:))));

% Quivers
doQ = 0;
nQx = 200;
nQz = 50;
[Z,X] = meshgrid(z,x);
ipxQ = fix(linspace(ipx1,ipx2,nQx));
ipzQ = fix(linspace(ipz1,ipz2,nQz));
[dataQx,dataQz] = meshgrid(ipxQ,ipzQ);
ipXQ = dataQx; ipZQ = dataQz;
dataQ.x = E.perp.x;
dataQ.y = E.perp.y;
dataQ.z = E.perp.z;
dataQ = J.perp;
maxQ = 2;
dataQ.abs = sqrt(dataQ.x.^2 + dataQ.z.^2);
dataQ.x(dataQ.abs>maxQ) = NaN;
dataQ.y(dataQ.abs>maxQ) = NaN;
dataQ.z(dataQ.abs>maxQ) = NaN;


% Panels
isub = 1;
tic;
for ivar = 1:nvars  
  hca = h(isub); isub = isub + 1;
  varstr = varstrs{ivar};
  variable = eval(varstr);  
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.Title.String = sprintf('%s',varstr); 
  hca.Title.Interpreter = 'none';  
  if any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  end
  hcb = colorbar('peer',hca);
  hb(isub-1) = hcb;
  %hcb.YLim = hca.CLim(2)*[-1 1];
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
%     for ixline = 1:size(saddle_locations,1)
%       sepA = saddle_values(ixline);
%       hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',sepA*[1 1],'color',cA.^4,'linewidth',2.0);  
%     end
    hold(hca,'off')  
  end
  if doQ
    hold(hca,'on')
    hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
    hold(hca,'off')  
  end
end


for ipanel = 1:npanels
  h(ipanel).YDir = 'normal';
  h(ipanel).XDir = 'reverse';
  h(ipanel).XLim = xlim;
  h(ipanel).YLim = zlim;
  if not(isempty(clim)), h(ipanel).CLim = clim; end
end
toc

%% Line plot at given x or z, define variable in cell array
% Define what variables to plot
varstrs = {'ve1.x','ve2.x','ve1.z','ve2.z','ve1.par','ve2.par','ExB.x','ExB.z','-ve1xB.x','-ve2xB.x','-ve1xB.z','-ve2xB.z','E.x','E.z'};
varstrs = {'ve1.x','ve2.x','B.z','E.z','-ve1xB.x','A'};
varstrs = {'ve1.perp.x','ve2.perp.x','vi1.perp.x','vi2.perp.x','ExB.x'};  
varstrs = {'ve1.perp.z','ve2.perp.z','vi1.perp.z','vi2.perp.z','ExB.z'};
varstrs = {'ne1','ne2','ni1','ni2','te2.scalar','ti2.scalar','pe2.scalar','pi2.scalar'};
varstrs = {'ne1','ne2','ni1','ni2'};
varstrs = {'Ute1','Ute2','Uti1','Uti2','Uke1','Uke2','Uki1','Uki2'}; clim = 12*[-1 1];
varstrs = {'pe1.xx','pe1.xy','pe1.yy','pe1.xz','pe1.zz','pe1.yz'}; clim = 0.25*[-1 1];
varstrs = {'ne1','ne2','ni1','ni2'}; clim = 2*[-1 1];
varstrs = {'pe2.xx','pe2.xy','pe2.yy','pe2.xz','pe2.zz','pe2.yz'}; clim = 0.25*[-1 1];
varstrs = {'pi2.xx','pi2.xy','pi2.yy','pi2.xz','pi2.zz','pi2.yz'}; clim = 0.3*[-1 1];
varstrs = {'ve1.x','ve1.y','ve1.z','ExB.x','ExB.y','ExB.z','-ve1xB.x','-ve1xB.y','-ve1xB.z','-E.x','-E.y','-E.z','-ve1xB.x-E.x','-ve1xB.y-E.y','-ve1xB.z-E.z'}; clim = [-1 1];
varstrs = {'vi1.x','vi1.y','vi1.z','ExB.x','ExB.y','ExB.z','vi1xB.x','vi1xB.y','vi1xB.z','E.x','E.y','E.z','vi1xB.x+E.x','vi1xB.y+E.y','vi1xB.z+E.z'}; clim = [-1 1];
varstrs = {'pi1.xx','pi1.xy','pi1.yy','pi1.xz','pi1.zz','pi1.yz','gradpi1.x','gradpi1.y','gradpi1.z'}; clim = 0.3*[-1 1];
varstrs = {'pi1.xx','pi1.zz','gradpi1.x','gradpi1.z','gradpi1_smooth.x','gradpi1_smooth.z','gradpi1.abs','gradpi1_smooth.abs'}; clim = 0.3*[-1 1];
varstrs = {'-ne1.*ve1xB.x','-ne1.*ve1xB.y','-ne1.*ve1xB.z','-ne1.*E.x','-ne1.*E.y','-ne1.*E.z','-ne1.*(ve1xB.x+E.x)','-ne1.*(ve1xB.y+E.y)','-ne1.*(ve1xB.z+E.z)','-gradpe1_smooth.x','-gradpe1_smooth.y','-gradpe1_smooth.z'}; clim = 0.5*[-1 1];
varstrs = {'-ne2.*ve2xB.x','-ne2.*ve2xB.y','-ne2.*ve2xB.z','-ne2.*E.x','-ne2.*E.y','-ne2.*E.z',...
           '-ne2.*(ve2xB.x+E.x)','-ne2.*(ve2xB.y+E.y)','-ne2.*(ve2xB.z+E.z)',...
           '-gradpe2_smooth.x','-gradpe2_smooth.y','-gradpe2_smooth.z',...
           '-ne2.*(ve2xB.x+E.x)-gradpe2_smooth.x','-ne2.*(ve2xB.y+E.y)-gradpe2_smooth.y','-ne2.*(ve2xB.z+E.z)-gradpe2_smooth.z'...
           }; clim = 0.5*[-1 1];
varstrs = {'ve1.x','ve2.x','vi1.x','vi2.x'}; clim = 3*[-1 1];
varstrs = {'te1.scalar','te2.scalar','te2.scalar./te1.scalar','ti1.scalar','ti2.scalar','ti2.scalar./ti1.scalar'}; clim = 0.8*[-1 1];
varstrs = {'te1.scalar','te2.scalar','te1.scalar./te2.scalar','ne1','ne2'}; clim = 3*[-1 1];
varstrs = {'ve1.par','ve1.perp.x','ve1.perp.y','ve1.perp.z','ve1.par./sqrt(ve1.perp.x.^2+ve1.perp.y.^2+ve1.perp.z.^2)'}; clim = 3*[-1 1];
%varstrs = {'vi1.par','vi1.perp.x','vi1.perp.y','vi1.perp.z','vi1.par./sqrt(vi1.perp.x.^2+vi1.perp.y.^2+vi1.perp.z.^2)'}; clim = 3*[-1 1];
varstrs = {'ne1','ne2'}; clim = 3*[-1 1];
varstrs = {'(-ve1xB.y+ve2xB.y)','(-vi1xB.y+vi2xB.y)','vi1.y','vi2.y'}; clim = 0.1*[-1 1];
varstrs = {'(-vi1xB.y+vi2xB.y)','vi1.y','vi2.y'}; clim = [];0.2*[-1 1];
varstrs = {'vi2xB.y','vi2xB.y_zx','vi2xB.y_xz'}; clim = [];0.2*[-1 1];
varstrs = {'vi1xB.y','vi2xB.y','vi1xB.y_zx','vi2xB.y_zx','vi1xB.y_xz','vi2xB.y_xz'}; clim = 0.5*[-1 1];
varstrs = {'vi1xB.y','vi2xB.y','-vi1xB.y+vi2xB.y','vi1xB.y_zx','vi2xB.y_zx','-vi1xB.y_zx+vi2xB.y_zx','vi1xB.y_xz','vi2xB.y_xz','-vi1xB.y_xz+vi2xB.y_xz'}; clim = 0.2*[-1 1];
%varstrs = {'-vi1xB.y+vi2xB.y','-vi1xB.y_zx+vi2xB.y_zx','-vi1xB.y_xz+vi2xB.y_xz','vi1.x','vi1.y','vi1.z','vi2.x','vi2.y','vi2.z'}; clim = 0.2*[-1 1];
%varstrs = {'vte1','vte2','vti1','vti2'}; clim = [];0.2*[-1 1];
%varstrs = {'wce1','wce2','wci1','wci2','vte1','vte2','vti1','vti2','re1','re2','ri1','ri2'}; clim = 10*[-1 1];0.2*[-1 1];
%varstrs = {'vte1','vte2','vti1','vti2','re1','re2','ri1','ri2'}; clim = 3*[-1 1];0.2*[-1 1];
varstrs = {'pi1.scalar','gradpi1.x','gradpi1.z','vExB.x','vExB.y','vExB.z'}; clim = [];
varstrs = {'vExB.x_yz','vExB.x_zy','vExB.y_zx','vExB.y_xz','vExB.z_xy','vExB.z_yx'}; clim = [];
varstrs = {'vExB.x','vExB.y','vExB.z','vExB.x_yz','vExB.x_zy','vExB.y_zx','vExB.y_xz','vExB.z_xy','vExB.z_yx'}; clim = [-2 2];
varstrs = {'vExB.x','vExB.x_yz','vExB.x_zy','vExB.y','vExB.y_zx','vExB.y_xz','vExB.z','vExB.z_xy','vExB.z_yx'}; clim = [-2 2];

varstrs = {'vExB.x','vExB.y','vExB.z','vDi1.x','vDi1.y','vDi1.z'}; clim = [];
varstrs = {'vExB.x','vExB.y','vExB.z','vDi1.x','vDi1.y','vDi1.z','vDi2.x','vDi2.y','vDi2.z','vDe1.x','vDe1.y','vDe1.z','vDe2.x','vDe2.y','vDe2.z'}; clim = [-1.1 1.1];
varstrs = {'vi1.y','vDi1.y','vExB.y','vDi1.y+vExB.y','vi2.y','vDi2.y','vExB.y','vDi2.y+vExB.y','ve1.y','vDe1.y','vExB.y','vDe1.y+vExB.y','ve2.y','vDe2.y','vExB.y','vDe2.y+vExB.y'}; clim = [-1.1 1.1];

varstrs = {{'B.z'},....
           {'ne1','ne2','ni1','ni2','ni1+ni2'},...
           {'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar','pB'},...
           {'te1.scalar','te2.scalar','ti1.scalar','ti2.scalar'},...
           {'ti1.scalar','ti2.scalar'}}; vallim = [];

varstrs = {{'B.z'},...
           {'ne1','ne2','ni1','ni2','ni1+ni2'},...
           ...%{'-(pe1.xx+pe1.xy+pe1.xz)','-(pe2.xx+pe2.xy+pe2.xz)','-(pi1.xx+pi1.xy+pi1.xz)','-(pi2.xx+pi2.xy+pi2.xz)'},...
           {'pe1.xx','pe1.xy','pe1.xz'},...
           {'pe2.xx','pe2.xy','pe2.xz'},...%{'0.5*(vve1.xx+vve1.xy+vve1.xz)','0.5*(vve2.xx+vve2.xy+vve2.xz)','0.5*(vvi1.xx+vvi1.xy+vvi1.xz)','0.5*(vvi2.xx+vvi2.xy+vvi2.xz)'},...
           {'pi1.xx','pi1.xy','pi1.xz'},...
           {'pi2.xx','pi2.xy','pi2.xz'}}; vallim = [];
if 1
  varstrs = {{'B.z'},...
             {'ne1','ne2','ni1','ni2','ni1+0.5*ni2'},...
             {'ve1.x','ve2.x','vi1.x','vi2.x'},...
             ...%{'0.5*(vve1.xx+vve1.xy+vve1.xz)','0.5*(vve2.xx+vve2.xy+vve2.xz)','0.5*(vvi1.xx+vvi1.xy+vvi1.xz)','0.5*(vvi2.xx+vvi2.xy+vvi2.xz)'},...
             ...%{'-0.5*vve1.xx','-0.5*vve2.xx','-0.5*vvi1.xx','-0.5*vvi2.xx'},...
             {'-pe1.xx','-pe2.xx','-pi1.xx','-pi2.xx'},...
             ...%{'pe1.xx+pe1.xy+pe1.xz','pe2.xx+pe2.xy+pe2.xz','pi1.xx+pi1.xy+pi1.xz','pi2.xx+pi2.xy+pi2.xz'},...
             ...%{'-pB','-0.5*BB.xx'},...%,'-0.5*BB.xy','-0.5*BB.xz'},...
             {'-0.5*B.abs.^2','-0.5*BB.xx'},...%,'-0.5*BB.xy','-0.5*BB.xz'},...
             ...%{'-pB','-0.5*(BB.xx+BB.xy+BB.xz)'},...
             ...%{'-0.5*(nmvve1.xx+nmvve1.xy+nmvve1.xz)','-0.5*(nmvve2.xx+nmvve2.xy+nmvve2.xz)','-0.5*(nmvvi1.xx+nmvvi1.xy+nmvvi1.xz)','-0.5*(nmvvi2.xx+nmvvi2.xy+nmvvi2.xz)'},...
             {'-0.5*nmvve1.xx','-0.5*nmvve2.xx','-0.5*nmvvi1.xx','-0.5*nmvvi2.xx'}}; vallim = [];
  colororders = {'1','bdacm','bdac','bdac','cd','bdac','','',''};
end
if 0
  varstrs = {{'B.x','E.z'},...
             {'ne1','ne2','ni1','ni2'},...
             {'ve1.x','ve2.x','vi1.x','vi2.x'},...
             {'ve1.y','ve2.y','vi1.y','vi2.y'},...
             {'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar'},...
             ...%{'pe1.xx+pe1.xy+pe1.xz','pe2.xx+pe2.xy+pe2.xz','pi1.xx+pi1.xy+pi1.xz','pi2.xx+pi2.xy+pi2.xz'},...
             ...%{'-pB','-0.5*(BB.xx+BB.xy+BB.xz)'},...
             {'0.5*B.abs.^2','pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar','0.5*B.abs.^2+pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar'}
             }; vallim = [];
  colororders = {'1y','bdacm','bdac','bdac','bdac','1ao','bdac','','',''};
end
if 0 % equilibrium
  varstrs = {{'B.x'},...
             {'ne1','ne2','ni1','ni2'},...
             {'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar'},...                          
             {'0.5*B.abs.^2','pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar','0.5*B.abs.^2+pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar'}
             }; vallim = [];
  colororders = {'1y','bdacm','bdac','1ao','bdac','','',''};
end
doAddSpeciesExplanation = 1;
speciesIdentification = {'i1','e1','i2','e2'};
speciesExplanation = {'hot ions','hot electrons','cold ions','cold electrons'};

doDiff = [];%[5 6 7];
var_operation = 'diff';
%doDiff = 1;

npanels = numel(varstrs);
nvars = cellfun(@numel, varstrs);

plotaxis = 'x'; % 'x' for horizontal cut, 'z' for vertical cut
zpick = 0;
xpick = 200;
zind = find_closest_ind(z,zpick);
xind = find_closest_ind(x,xpick);

%xlim = x([1 end/2])'+[100 -100];
xlim = [150 x(fix(end/2))];
xlim = [x(fix(end/2)) x(end)-150];
xlim = x([1 end])'+[150 -150];
zlim = [-10 10];

ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
%ipz1:2:ipz2;
%[X,Z] = ndgrid(x,z);
%Xp = X(ipx,ipz);
%Zp = Z(ipx,ipz);
switch plotaxis
  case 'x'
    ipx = ipx1:2:ipx2;
    ipz = zind;
    plot_dep = x(ipx);
    plotlim = xlim;
    pickind = zpick;
    pickval = z(ipz);
    pickstr = 'z';
  case 'z'
    ipz = ipz1:2:ipz2;
    ipx = xind;
    plot_dep = z(ipz);
    plotlim = zlim;
    pickind = xpick;
    pickval = x(ipx);
    pickstr = 'x';
end

linewidth = 1.5;
fontsize = 12;

% Initialize figure
npanels = npanels;
nrows = npanels;
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
isub = 1; 
clear h hleg;
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1; 
  hold(h(isub-1),'off')

  %h(isub).Box = 'on';
end

    
% Panels
isub = 1;
tic;
doColor = 0;
for ipanel = 1:npanels
  hca = h(isub); isub = isub + 1;
  colors = pic_colors(colororders{ipanel});
  if size(colors,1)>=(nvars(ipanel)), doColor = 1; else doColor = 0; end
  hplot = [];
  for ivar = 1:nvars(ipanel)  
    if ivar == 1, hold(hca,'off');
    elseif ivar == 2,  hold(hca,'on'); end 
    varstr = varstrs{ipanel}{ivar};
    variable = eval(varstr);
    
    if intersect(doDiff,ipanel)
      diff_variable = [0; diff(variable(ipx,ipz))];
      hplot_tmp = plot(hca,plot_dep,diff_variable,'LineWidth',linewidth);
      varstrs{ipanel}{ivar} = sprintf('diff(%s)',varstrs{ipanel}{ivar});      
    else
      hplot_tmp = plot(hca,plot_dep,variable(ipx,ipz),'LineWidth',linewidth);      
    end
    if doColor, hplot_tmp.Color = colors(ivar,:); end
    hplot{ivar} = hplot_tmp;
    hca.XLabel.String = sprintf('%s (d_i)',plotaxis);
    %hca.YLabel.String = varstr;
    %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
    %hca.Title.String = sprintf('%s',varstr); 
    %hca.YLabel.Interpreter = 'none';        
  end
  leg_str_tmp = varstrs{ipanel};
  if doAddSpeciesExplanation
    for ivar = 1:nvars(ipanel)  
      for ispecies = 1:numel(speciesIdentification)
        if strfind(leg_str_tmp{ivar},speciesIdentification{ispecies})
          leg_str_tmp{ivar} = sprintf('%s (%s)',leg_str_tmp{ivar},speciesExplanation{ispecies});
        end
      end    
    end
  end
  hleg_tmp = legend(hca,leg_str_tmp,'box','off','fontsize',fontsize,'location','eastoutside');
  hleg(ipanel) = hleg_tmp;
  hca.FontSize = fontsize;
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XDir = 'reverse';
end
%hlink = linkprop(h(4:6),'YLim');
hlink = linkprop(h,'XLim'); 
if strcmp(plotaxis,'x'), h(1).XLim = xlim;
elseif strcmp(plotaxis,'z'), h(1).XLim = zlim;
end
%hlink = linkprop(hleg,'Position');
%addprop(hlink,'PlotBoxAspectRatio')
h(1).Title.String = sprintf('%s = %.0f (d_i), time = %g (1/wci) = %g (1/wpe)',pickstr,pickval,time,timestep);
arrayfun(@(x)eval(sprintf('x.Position(3) = 0.6;'),x),h)
drawnow
compact_panels(0.012)

%% (OLD) Plot, 4 species plasma properties, 1 species per column
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