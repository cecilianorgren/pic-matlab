function out = function_load_data(txtfile)
% FUNCTION_LOAD_DATA Loads data from particle-in-cell simulation.
%   Loads data from .dat file, which were created using Fortran 90.   
%
% Example:
%   filepath = '/Path/to/your/file/fields-24000.dat';
%   data = function_load_data(filepath);
%
%   hca = subplot(3,1,1);
%   pcolor(hca,data.x_di,data.z_di,data.Ey')
%   shading(hca,'flat')
%   hca.XLabel.String = 'x (d_i)';
%   hca.YLabel.String = 'z (d_i)';
%   hcb = colorbar('peer',hca);
%   hcb.YLabel.String = 'E_y (v_AB_0)';
%
%   hca = subplot(3,1,2);
%   pcolor(hca,data.x_di,data.z_di,data.cold_ion_top.n')
%   shading(hca,'flat')
%   hca.XLabel.String = 'x (d_i)';
%   hca.YLabel.String = 'z (d_i)';
%   hcb = colorbar('peer',hca);
%   hcb.YLabel.String = 'n_{i,cold,top} (n_0)';
%
%   hca = subplot(3,1,3);
%   pcolor(hca,data.x_di,data.z_di,data.cold_ion.n')
%   shading(hca,'flat')
%   hca.XLabel.String = 'x (d_i)';
%   hca.YLabel.String = 'z (d_i)';
%   hcb = colorbar('peer',hca);
%   hcb.YLabel.String = 'n_{i,cold} (n_0)';

%% Read data from fields-_____.dat file

nss = 6; % Oxygen run has 4 species
[fid, message] = fopen(txtfile,'r','ieee-le');
      if fid < 0
        error('Failed to open file "%s" because "%s"', txtfile, message);
      end
header = fread(fid,1,'integer*8');

it = fread(fid,1,'integer*4');                                % it
dt = fread(fid,1,'real*4');                                   % dt
teti = fread(fid,1,'real*4');                                 % teti
xmax = fread(fid,1,'real*4');                                 % xmax
zmax = fread(fid,1,'real*4');                                 % zmax
nnx = fread(fid,1,'integer*4');                               % nnx
nnz = fread(fid,1,'integer*4');                               % nnz

vxs = zeros(nnx,nnz,nss);
vys = zeros(nnx,nnz,nss);
vzs = zeros(nnx,nnz,nss);

for is = 1:nss, vxs(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxs 
for is = 1:nss, vys(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vys 
for is = 1:nss, vzs(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vzs
bx = fread(fid,[nnx nnz],'real*4');                           % bx 
by = fread(fid,[nnx nnz],'real*4');                           % by 
bz = fread(fid,[nnx nnz],'real*4');                           % bz 
ex = fread(fid,[nnx nnz],'real*4');                           % ex 
ey = fread(fid,[nnx nnz],'real*4');                           % ey 
ez = fread(fid,[nnx nnz],'real*4');                           % ez 

dns = zeros(nnx,nnz,nss);
for is = 1:nss, dns(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % dns 

xe = fread(fid,nnx,'real*4');                                 % xe 
ze = fread(fid,nnz,'real*4');                                 % ze 
mass = fread(fid,nss,'real*4');                               % mass 
q = fread(fid,nss,'real*4');                                  % q 
time = fread(fid,1,'real*8');                                 % time 
wpewce = fread(fid,1,'real*4');                               % wpewce 
dfac = fread(fid,nss,'real*4');                               % dfac

vxx = zeros(nnx,nnz,nss);
vyy = zeros(nnx,nnz,nss);
vzz = zeros(nnx,nnz,nss);
vxy = zeros(nnx,nnz,nss);
vxz = zeros(nnx,nnz,nss);
vyz = zeros(nnx,nnz,nss);

for is = 1:nss, vxx(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxx 
for is = 1:nss, vyy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyy 
for is = 1:nss, vzz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pzz 
for is = 1:nss, vxy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxy 
for is = 1:nss, vxz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxz 
for is = 1:nss, vyz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyz
remainder = fread(fid); 

%% Normalize data
mime = mass(1)/mass(2);
ex = ex*sqrt(mime)*wpewce^2;
ey = ey*sqrt(mime)*wpewce^2;
ez = ez*sqrt(mime)*wpewce^2;
bx = bx*wpewce;
by = by*wpewce;
bz = bz*wpewce;

iSpecies = 1; % Hot ions
data_hot_ions           = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = 2; % Hot electrons
data_hot_electrons      = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = 3; % Cold ions, top
data_cold_ions_top      = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = 4; % Cold electrons, top
data_cold_electrons_top = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = 5; % Cold ions, bottom
data_cold_ions_bot      = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = 6; % Cold electrons, bottom
data_cold_electrons_bot = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = [3 5]; % Cold ions, top and bottom
data_cold_ions          = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

iSpecies = [4 6]; % Cold electrons, top and bottom
data_cold_electrons     = load_species(iSpecies,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz);

% Collect data into structure
out.species = {'hot ions','hot electrons','cold ions top','cold electrons top','cold ions bot','cold electrons bot'};
out.time_wpe = time; % time units of inverse electron plasma frequencies
out.time_wci = time/(wpewce*mime); % time units of inverse ion cyclotron frequencies
out.nx = nnx;
out.nz = nnz;
out.x_de = xe;
out.z_de = ze;
out.x_di = xe/sqrt(mime);
out.z_di = ze/sqrt(mime);
out.mass = mass;
out.q = q;
out.wpewce = wpewce;
out.Ex = ex;
out.Ey = ey;
out.Ez = ez;
out.Bx = bx;
out.By = by;
out.Bz = bz;
out.hot_ion = data_hot_ions;
out.hot_ele = data_hot_ions;
out.cold_ion_top = data_cold_ions_top;
out.cold_ion_bot = data_cold_ions_bot;
out.cold_ele_top = data_cold_electrons_top;
out.cold_ele_bot = data_cold_electrons_bot;
out.cold_ion = data_cold_ions;
out.cold_ele = data_cold_electrons;

function out = load_species(species,nnx,nnz,mass,dfac,wpewce,dns,vxs,vys,vzs,vxx,vxy,vxz,vyy,vyz,vzz)
  % LOAD_SPECIES Normalizes data for given species.  
  mime = mass(1)/mass(2);
  mass = mass(species(1));
  
  % Initialize variables
  n_tmp = zeros(nnx,nnz);
  vxs_tmp = zeros(nnx,nnz);
  vys_tmp = zeros(nnx,nnz);
  vzs_tmp = zeros(nnx,nnz);
  vxx_tmp = zeros(nnx,nnz);
  vxy_tmp = zeros(nnx,nnz);
  vxz_tmp = zeros(nnx,nnz);
  vyy_tmp = zeros(nnx,nnz);
  vyz_tmp = zeros(nnx,nnz);
  vzz_tmp = zeros(nnx,nnz);
          
  % Collect data for given species
  for iSpecies = species
    % density      
   n_tmp = n_tmp + dns(:,:,iSpecies)*dfac(iSpecies);
    % flux
    vxs_tmp = vxs_tmp + vxs(:,:,iSpecies)*dfac(iSpecies);
    vys_tmp = vys_tmp + vys(:,:,iSpecies)*dfac(iSpecies);
    vzs_tmp = vzs_tmp + vzs(:,:,iSpecies)*dfac(iSpecies);
     % stress tensor
    vxx_tmp = vxx_tmp + vxx(:,:,iSpecies)*dfac(iSpecies);
    vxy_tmp = vxy_tmp + vxy(:,:,iSpecies)*dfac(iSpecies);
    vxz_tmp = vxz_tmp + vxz(:,:,iSpecies)*dfac(iSpecies);
    vyy_tmp = vyy_tmp + vyy(:,:,iSpecies)*dfac(iSpecies);
    vyz_tmp = vyz_tmp + vyz(:,:,iSpecies)*dfac(iSpecies);
    vzz_tmp = vzz_tmp + vzz(:,:,iSpecies)*dfac(iSpecies);
  end        
  
  % Normalize and collect data into a structure
  % density
  out.n = n_tmp;
  % flux 
  out.jx = vxs_tmp*wpewce*sqrt(mime);
  out.jy = vys_tmp*wpewce*sqrt(mime);
  out.jz = vzs_tmp*wpewce*sqrt(mime);
  % velocity 
  out.vx = out.jx./out.n;
  out.vy = out.jy./out.n;
  out.vz = out.jz./out.n;
  % pressure
  out.pxx = (vxx_tmp - vxs_tmp.*vxs_tmp./n_tmp)*mass*wpewce^2;
  out.pxy = (vxy_tmp - vxs_tmp.*vys_tmp./n_tmp)*mass*wpewce^2;
  out.pxz = (vxz_tmp - vxs_tmp.*vzs_tmp./n_tmp)*mass*wpewce^2;
  out.pyy = (vyy_tmp - vys_tmp.*vys_tmp./n_tmp)*mass*wpewce^2;
  out.pyz = (vyz_tmp - vys_tmp.*vzs_tmp./n_tmp)*mass*wpewce^2;
  out.pzz = (vzz_tmp - vzs_tmp.*vzs_tmp./n_tmp)*mass*wpewce^2;
  out.p = (out.pxx + out.pyy + out.pzz)/3;
  % temperature
  out.t = out.p./out.n;
  
