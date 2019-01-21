% load file
nsp = 4; % Run has 4 species
listFiles = dir(pathParticles);
nf = 1;
nameFile = sprintf('%.0f.dat',nf);

nss = 4; % number of species
numberel = 101; 

% preallocate
fxy = zeros(numberel,numberel,nss);
fxz = zeros(numberel,numberel,nss);
fyz = zeros(numberel,numberel,nss);
axes = zeros(numberel,nss);

fid = fopen([pathParticles nameFile],'r','ieee-le');
header = fread(fid,1,'integer*8');
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
 
fclose(fid);
%remainder = fread(fid);                                        % what's left of the file