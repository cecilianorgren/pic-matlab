function [xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
    jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz, ...
    wpewce,mass,vex,vey,vez, pxxi,pyyi,pzzi,pxxe,pyye,pzze,pi,pe,pyzi,pxyi,a,pxye,pyze] ...
    = func_read_file_fields_invar(txtfile);

%% Read data from fields-_____.dat file
tic
nss = 4; % Oxygen run has 4 species
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

% bx = zeros(nnx,nnz);
% by = zeros(nnx,nnz);
% bz = zeros(nnx,nnz);
% ex = zeros(nnx,nnz);
% ey = zeros(nnx,nnz);
% ez = zeros(nnx,nnz);

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
pxx = zeros(nnx,nnz,nss);
pyy = zeros(nnx,nnz,nss);
pzz = zeros(nnx,nnz,nss);
pxy = zeros(nnx,nnz,nss);
pxz = zeros(nnx,nnz,nss);
pyz = zeros(nnx,nnz,nss);

for is = 1:nss, pxx(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxx 
for is = 1:nss, pyy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyy 
for is = 1:nss, pzz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pzz 
for is = 1:nss, pxy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxy 
for is = 1:nss, pxz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxz 
for is = 1:nss, pyz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyz
remainder = fread(fid); 
toc


%% Finishing reading fields file
threshold=0;


dni=0;dne=0;jix=0;jiy=0;jiy=0;jiz=0;jex=0;jey=0;jez=0;pxxi=0;pyyi=0;pzzi=0;
pxyi=0;pxzi=0;pxxe=0;pyye=0;pzze=0;pxye=0;pxze=0;pyzi=0;pyze=0;
%% Normalize lengths and fields to V_A and B_0 and c/omega_pi

xe=xe/sqrt(mass(1));
ze=ze/sqrt(mass(1));
bx=bx*wpewce(1);
by=by*wpewce(1);
bz=bz*wpewce(1);
ex=ex*sqrt(mass(1))*wpewce(1)^2;
ey=ey*sqrt(mass(1))*wpewce(1)^2;
ez=ez*sqrt(mass(1))*wpewce(1)^2;

for isp=1:2:nss
    ion = isp;
    ele = isp+1;
    
    dni=dni+squeeze(dns(:,:,ion(1)))*dfac(ion(1));
    dne=dne+squeeze(dns(:,:,ele(1)))*dfac(ele(1));
    
    jix=jix+squeeze(vxs(:,:,ion(1)))*dfac(ion(1));
    jiy=jiy+squeeze(vys(:,:,ion(1)))*dfac(ion(1));
    jiz=jiz+squeeze(vzs(:,:,ion(1)))*dfac(ion(1));    
    jex=jex+squeeze(vxs(:,:,ele(1)))*dfac(ele(1)) ;
    jey=jey+squeeze(vys(:,:,ele(1)))*dfac(ele(1)) ;
    jez=jez+squeeze(vzs(:,:,ele(1)))*dfac(ele(1));
    
    pxxi=pxxi+squeeze(pxx(:,:,ion))*dfac(ion);
    pyyi=pyyi+squeeze(pyy(:,:,ion))*dfac(ion);
    pzzi=pzzi+squeeze(pzz(:,:,ion))*dfac(ion);
    pxyi=pxyi+squeeze(pxy(:,:,ion))*dfac(ion);
    pxzi=pxzi+squeeze(pxz(:,:,ion))*dfac(ion);
    pyzi=pyzi+squeeze(pyz(:,:,ion))*dfac(ion);
    
    pxxe=pxxe+squeeze(pxx(:,:,ele))*dfac(ele);
    pyye=pyye+squeeze(pyy(:,:,ele))*dfac(ele);
    pzze=pzze+squeeze(pzz(:,:,ele))*dfac(ele);
    pxye=pxye+squeeze(pxy(:,:,ele))*dfac(ele);
    pxze=pxze+squeeze(pxz(:,:,ele))*dfac(ele);
    pyze=pyze+squeeze(pyz(:,:,ele))*dfac(ele);    
    
end

% Set species back to primary indicies 
ion = 1; ele = 2;


tmpvix = jix./dni; tmpvix(dni < threshold) = 0;
tmpviy = jiy./dni; tmpviy(dni < threshold) = 0;
tmpviz = jiz./dni; tmpviz(dni < threshold) = 0;

    pxzi=mass(ion)*(pxzi-tmpvix.*jiz);
    pxxi=mass(ion)*(pxxi-tmpvix.*jix);
    pxyi=mass(ion)*(pxyi-tmpvix.*jiy);
    pyyi=mass(ion)*(pyyi-tmpviy.*jiy);
    pyzi=mass(ion)*(pyzi-tmpviy.*jiz);
    pzzi=mass(ion)*(pzzi-tmpviz.*jiz);

tmpvex = jex./dne; tmpvex(dne < threshold) = 0;
tmpvey = jey./dne; tmpvey(dne < threshold) = 0;
tmpvez = jez./dne; tmpvez(dne < threshold) = 0;    
    
    pxze=mass(ele)*(pxze-tmpvex.*jez);
    pxxe=mass(ele)*(pxxe-tmpvex.*jex);
    pxye=mass(ele)*(pxye-tmpvex.*jey);
    pyye=mass(ele)*(pyye-tmpvey.*jey);
    pyze=mass(ele)*(pyze-tmpvey.*jez);
    pzze=mass(ele)*(pzze-tmpvez.*jez);
    

%% Normalize derived quantities

    pxxi=pxxi*wpewce(1)^2;
    pyyi=pyyi*wpewce(1)^2;
    pzzi=pzzi*wpewce(1)^2;
    pxyi=pxyi*wpewce(1)^2;
    pxzi=pxzi*wpewce(1)^2;
    pyzi=pyzi*wpewce(1)^2;
    
    pxxe=pxxe*wpewce(1)^2;
    pyye=pyye*wpewce(1)^2;
    pzze=pzze*wpewce(1)^2;
    pxye=pxye*wpewce(1)^2;
    pxze=pxze*wpewce(1)^2;
    pyze=pyze*wpewce(1)^2;

    jix=jix*wpewce(1)*sqrt(mass(1));
    jiy=jiy*wpewce(1)*sqrt(mass(1));
    jiz=jiz*wpewce(1)*sqrt(mass(1));
    
    jex=jex*wpewce(1)*sqrt(mass(1));
    jey=jey*wpewce(1)*sqrt(mass(1));
    jez=jez*wpewce(1)*sqrt(mass(1));
    pi=(pxxi+pyyi+pzzi)/3;
    pe=(pxxe+pyye+pzze)/3;
    
    %% from these normalized quantities we can derive additional quantities
    ti = pi./dni; ti(dni<threshold)=0;
    te = pe./dne; te(dne<threshold)=0;
% michaels invar.pro stops here.... 


    vix = jix./dni; vix(dni < threshold) = 0;
    viy = jiy./dni; viy(dni < threshold) = 0;
    viz = jiz./dni; viz(dni < threshold) = 0;
    
    vex = jex./dne; vex(dne < threshold) = 0;
    vey = jey./dne; vey(dne < threshold) = 0;
    vez = jez./dne; vez(dne < threshold) = 0; 
    
 dz = ze(2)-ze(1);
 dx = xe(2)-xe(1);
ixm = 10;
a = zeros(nnx,nnz);
a(ixm,1) = 0;
      for j=2:nnz 
        a(ixm,j) = a(ixm,j-1) + dz*bx(ixm,j-1);
      end
      for ix=ixm+1:nnx
        a(ix,:) = a(ix-1,:) - bz(ix-1,:)*dx; % ;     advance to the right
      end
      for ix=ixm-1:-1:1
        a(ix,:) = a(ix+1,:)+ bz(ix,:)*dx; % ;     advance to the left
      end    
    
    
    
    
