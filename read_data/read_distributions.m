%% Read distribution function binary files
 
function [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
    = read_distributions(txtfile,nss,nv)

fid = fopen(txtfile,'r','ieee-le');

%nss = 4; % Oxygen run has 4 species

numberel = nv; %hardcoded? 

%% preallocate
fxy = zeros(numberel,numberel,nss);
fxz = zeros(numberel,numberel,nss);
fyz = zeros(numberel,numberel,nss);
axes = zeros(numberel,nss);

%% read file
fread(fid,1,'integer*8');
for is = 1:nss, axes(:,is) = fread(fid,numberel,'real*4'); end
xlo = fread(fid,1,'real*4');                                 
xhi = fread(fid,1,'real*4');                                 
zlo = fread(fid,1,'real*4');
zhi = fread(fid,1,'real*4');
ic = fread(fid,nss,'integer*4'); % long
fxyz=fread(fid,numberel*numberel*numberel*nss,'real*4');
fxyz=reshape(fxyz,numberel,numberel,numberel,nss);
for is = 1:nss, fxy(:,:,is) = fread(fid,[numberel numberel],'real*4'); end  % vxs 
for is = 1:nss, fxz(:,:,is) = fread(fid,[numberel numberel],'real*4'); end  % vxs 
for is = 1:nss, fyz(:,:,is) = fread(fid,[numberel numberel],'real*4'); end  % vxs 
vxa = fread(fid,nss,'real*4');
vya = fread(fid,nss,'real*4');
vza = fread(fid,nss,'real*4');

remainder = fread(fid);
%numel(remainder)

%% readu,1,axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza

% AXES (VERT)     FLOAT     = Array[101, 4]
% FXY (VERT)      FLOAT     = Array[101, 101, 4]
% FXYZ (VERT)     FLOAT     = Array[101, 101, 101, 4]
% FXZ (VERT)      FLOAT     = Array[101, 101, 4]
% FYZ (VERT)      FLOAT     = Array[101, 101, 4]
% IC (VERT)       LONG      = Array[4]
% INDISTF         UNDEFINED = <Undefined>
% NDIST           INT       =      101
% NSS             INT       =        4
% NX              INT       =     3200
% NZ              INT       =     1600
% VXA (VERT)      FLOAT     = Array[4]
% VYA (VERT)      FLOAT     = Array[4]
% VZA (VERT)      FLOAT     = Array[4]
% XHI (VERT)      FLOAT     = Array[1]
% XLO (VERT)      FLOAT     = Array[1]
% ZHI (VERT)      FLOAT     = Array[1]
% ZLO (VERT)      FLOAT     = Array[1]
fclose(fid);