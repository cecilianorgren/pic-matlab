% What the output should be, used for finding header
% cd /Users/cno062/MATLAB/pic-matlab/
% system('./read_fields.exe')

% read(20)it,dt,teti,xmax,zmax,nnx,nnz,vxs,vys,vzs,
%      &         bx,by,bz,ex,ey,ez,dns,xe,ze,mass,q,time,wpewce,dfac,
%      &         pxx,pyy,pzz,pxy,pxz,pyz

% load file
nsp = 4; % Run has 4 species
fid = fopen([pathFields 'fields-05978.dat'],'r','ieee-le');

%% Read data from fields-_____.dat file
header = fread(fid,1,'integer*8');

it = fread(fid,1,'integer*4');                                % it
dt = fread(fid,1,'real*4');                                   % dt
teti = fread(fid,1,'real*4');                                 % teti
xmax = fread(fid,1,'real*4');                                 % xmax
zmax = fread(fid,1,'real*4');                                 % zmax
nnx = fread(fid,1,'integer*4');                               % nnx
nnz = fread(fid,1,'integer*4');                               % nnz
for is = 1:nsp, vxs(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxs 
for is = 1:nsp, vys(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vys 
for is = 1:nsp, vzs(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vzs
bx = fread(fid,[nnx nnz],'real*4');                           % bx 
by = fread(fid,[nnx nnz],'real*4');                           % by 
bz = fread(fid,[nnx nnz],'real*4');                           % bz 
ex = fread(fid,[nnx nnz],'real*4');                           % ex 
ey = fread(fid,[nnx nnz],'real*4');                           % ey 
ez = fread(fid,[nnx nnz],'real*4');                           % ez 
for is = 1:nsp, dns(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % dns 
xe = fread(fid,nnx,'real*4');                                 % xe 
ze = fread(fid,nnz,'real*4');                                 % ze 
mass = fread(fid,nsp,'real*4');                               % mass 
q = fread(fid,nsp,'real*4');                                  % q 
time = fread(fid,1,'real*8');                                 % time 
wpewce = fread(fid,1,'real*4');                               % wpewce 
dfac = fread(fid,nsp,'real*4');                               % dfac
for is = 1:nsp, pxx(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxx 
for is = 1:nsp, pyy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyy 
for is = 1:nsp, pzz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pzz 
for is = 1:nsp, pxy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxy 
for is = 1:nsp, pxz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxz 
for is = 1:nsp, pyz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyz
remainder = fread(fid) 

%% Read data from fields-_____.dat file (c_eval)
header = fread(fid,1,'integer*8');

it = fread(fid,1,'integer*4');                                % it
dt = fread(fid,1,'real*4');                                   % dt
teti = fread(fid,1,'real*4');                                 % teti
xmax = fread(fid,1,'real*4');                                 % xmax
zmax = fread(fid,1,'real*4');                                 % zmax
nnx = fread(fid,1,'integer*4');                               % nnx
nnz = fread(fid,1,'integer*4');                               % nnz
c_eval('vxs(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % vxs 
c_eval('vys(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % vys 
c_eval('vzs(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % vzs
bx = fread(fid,[nnx nnz],'real*4');                           % bx 
by = fread(fid,[nnx nnz],'real*4');                           % by 
bz = fread(fid,[nnx nnz],'real*4');                           % bz 
ex = fread(fid,[nnx nnz],'real*4');                           % ex 
ey = fread(fid,[nnx nnz],'real*4');                           % ey 
ez = fread(fid,[nnx nnz],'real*4');                           % ez 
c_eval('dns? = fread(fid,[nnx nnz],''real*4'');',1:nsp)       % dns 
xe = fread(fid,nnx,'real*4');                                 % xe 
ze = fread(fid,nnz,'real*4');                                 % ze 
mass = fread(fid,nsp,'real*4');                               % mass 
q = fread(fid,nsp,'real*4');                                  % q 
time = fread(fid,1,'real*8');                                 % time 
wpewce = fread(fid,1,'real*4');                               % wpewce 
dfac = fread(fid,nsp,'real*4');                               % dfac
c_eval('pxx(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % pxx 
c_eval('pyy(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % pyy 
c_eval('pzz(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % pzz 
c_eval('pxy(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % pxy 
c_eval('pxz(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % pxz 
c_eval('pyz(:,:,?) = fread(fid,[nnx nnz],''real*4'');',1:nsp) % pyz
remainder = fread(fid)                                        % what's left of the file


%% Plot data
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

isub = 1;
isp = 2;
if 1 % vxs
  hca = h(isub); isub = isub + 1;
  pcolor(hca,xe,ze,vxs(:,:,isp)'); 
  shading(hca,'flat')
  hca.Title.String = 'vxs2';
  colorbar('peer',hca)
end
if 1 % vys
  hca = h(isub); isub = isub + 1;
  pcolor(hca,xe,ze,vys(:,:,isp)'); 
  shading(hca,'flat')
  hca.Title.String = 'vys2';
  colorbar('peer',hca)
end
if 1 % vzs
  hca = h(isub); isub = isub + 1;
  pcolor(hca,xe,ze,vzs(:,:,isp)'); 
  shading(hca,'flat')
  hca.Title.String = 'vzs2';
  colorbar('peer',hca)
end
if 1 % bx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,bx'); 
  shading(hca,'flat')
  hca.Title.String = 'bx';
  colorbar('peer',hca)
end
if 1 % by
  hca = h(isub); isub = isub + 1;
  pcolor(hca,by'); 
  shading(hca,'flat')
  hca.Title.String = 'by';
  colorbar('peer',hca)
end
if 1 % bz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,bz'); 
  shading(hca,'flat')
  hca.Title.String = 'bz';
  colorbar('peer',hca)
end
if 1 % ex
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ex'); 
  shading(hca,'flat')
  hca.Title.String = 'ex';
  colorbar('peer',hca)
end
if 1 % ey
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ey'); 
  shading(hca,'flat')
  hca.Title.String = 'ey';
  colorbar('peer',hca)
end
if 1 % ez
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ez'); 
  shading(hca,'flat')
  hca.Title.String = 'ez';
  colorbar('peer',hca)
end

cmap = irf_colormap('poynting');
for ip = 1:npanels
  hca = h(ip);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  %colormap(hca,cn.cmap('bluered3'))
  colormap(hca,cmap)
end
