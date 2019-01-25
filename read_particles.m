function [xe,ze,e,b,ni1,ne1,ni2,ne2,vi1,ve1,vi2,ve2,ji1,je1,ji2,je2,...
    pi1,pe1,pi2,pe2,ti1,te1,ti2,te2,dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] ...
    = read_fields(txtfile,varargin)  
  % READ_FILES
  %   Reads files and return normalized quantities
 
  %% defaults
  nss = 4; % numer of species
  numberel = 101; % number of bins for particle distributions
  %groups = {[1 3],[2 4]}; % electron and ions
  
  % read input
  nargs = numel(varargin);
  have_options = nargs > 1;
  while have_options
    switch(lower(args{1}))
      case {'nss'} % number of species
        l = 2;
        nss = args{2};
      case 'numberel' % number of bins for particle distributions
        l = 2;
        numberel = args{2};
%       case 'groups'
%         l = 2;
%         groups = args{2};
%         doGroupManual = 1;
%       case 'groupcharge'
%         l = 1;
%         doGroupCharge = 1;
%       case 'groupmass'
%         l = 1;
%         doGroupMass = 1;
    end
    args = args((l+1):end);
    if isempty(args), break, end
  end
    
%   if doGroupManual == 1, doGroupCharge = 0; doGroupMass = 0;    
%   elseif doGroupMass == 1, doGroupCharge = 0;      
%   end
  
  %% load data
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
 
  
  % Set groups  
%   if doGroupMass
%     uniqueMass = unique(mass);
%     for imass = 1:numel(uniqueMass)
%       groups{imass} = find(mass==uniqueMAss(imass));
%     end
%   elseif doGroupCharge
%     uniqueCharge = unique(q);
%     for icharge = 1:numel(uniqueCharge)
%       groups{icharge} = find(q==uniqueCharge(icharge));
%     end
%   end
  
  %% Normalize 
  % 
  % readu,1,it,dt,teti,xmax,zmax,nx0,nz0,vxs,vys,vzs,    $
  %                bx,by,bz,ex,ey,ez,dns,x,z,mass,q,c,wpewce,dfac,$
  %               pxx,pyy,pzz,pxy,pxz,pyz

  xe=xe/sqrt(mass(1));
  ze=ze/sqrt(mass(1));
  bx=bx*wpewce(1);
  by=by*wpewce(1);
  bz=bz*wpewce(1);
  ex=ex*sqrt(mass(1))*wpewce(1)^2;
  ey=ey*sqrt(mass(1))*wpewce(1)^2;
  ez=ez*sqrt(mass(1))*wpewce(1)^2;
  
  time = time/wpewce/sqrt(mass(1));
  
%% Species 1,2
ion = 1;
ele = 2;

dni=squeeze(dns(:,:,ion(1)))*dfac(ion(1));
dne=squeeze(dns(:,:,ele(1)))*dfac(ele(1));

jix=squeeze(vxs(:,:,ion(1)))*dfac(ion(1));
jiy=squeeze(vys(:,:,ion(1)))*dfac(ion(1));
jiz=squeeze(vzs(:,:,ion(1)))*dfac(ion(1));



vix = jix./dni;
viy = jiy./dni;
viz = jiz./dni;
vix(dni < 0.005) = 0;
viy(dni < 0.005) = 0;
viz(dni < 0.005) = 0;

jex=squeeze(vxs(:,:,ele(1)))*dfac(ele(1)) ;
jey=squeeze(vys(:,:,ele(1)))*dfac(ele(1)) ;
jez=squeeze(vzs(:,:,ele(1)))*dfac(ele(1));

vex = jex./dne;
vey = jey./dne;
vez = jez./dne;

vex(dne < 0.005) = 0;
vey(dne < 0.005) = 0;
vez(dne < 0.005) = 0;

pxxi=squeeze(pxx(:,:,ion))*dfac(ion);
pyyi=squeeze(pyy(:,:,ion))*dfac(ion);
pzzi=squeeze(pzz(:,:,ion))*dfac(ion);
pxyi=squeeze(pxy(:,:,ion))*dfac(ion);
pxzi=squeeze(pxz(:,:,ion))*dfac(ion);
pyzi=squeeze(pyz(:,:,ion))*dfac(ion);

pxzi=mass(ion)*(pxzi-vix.*jiz);
pxxi=mass(ion)*(pxxi-vix.*jix);
pxyi=mass(ion)*(pxyi-vix.*jiy);
pyyi=mass(ion)*(pyyi-viy.*jiy);
pyzi=mass(ion)*(pyzi-viy.*jiz);
pzzi=mass(ion)*(pzzi-viz.*jiz);

% Move to right frame

%     pxzi=mass(ion)*(pxzi-vix.*jiz-jix.*viz+dni.*vix.*viz);
%     pxxi=mass(ion)*(pxxi-vix.*jix-jix.*vix+dni.*vix.*vix);
%     pxyi=mass(ion)*(pxyi-vix.*jiy-jix.*viy+dni.*vix.*viy);
%     pyyi=mass(ion)*(pyyi-viy.*jiy-jiy.*viy+dni.*viy.*viy);
%     pyzi=mass(ion)*(pyzi-viy.*jiz-jiy.*viz+dni.*viy.*viz);
%     pzzi=mass(ion)*(pzzi-viz.*jiz-jiz.*viz+dni.*viz.*viz);


pxxe=squeeze(pxx(:,:,ele))*dfac(ele);
pyye=squeeze(pyy(:,:,ele))*dfac(ele);
pzze=squeeze(pzz(:,:,ele))*dfac(ele);
pxye=squeeze(pxy(:,:,ele))*dfac(ele);
pxze=squeeze(pxz(:,:,ele))*dfac(ele);
pyze=squeeze(pyz(:,:,ele))*dfac(ele);

pxze=mass(ele)*(pxze-vex.*jez);
pxxe=mass(ele)*(pxxe-vex.*jex);
pxye=mass(ele)*(pxye-vex.*jey);
pyye=mass(ele)*(pyye-vey.*jey);
pyze=mass(ele)*(pyze-vey.*jez);
pzze=mass(ele)*(pzze-vez.*jez);

%     pxze=mass(ele)*(pxze-vex.*jez-jex.*vez+dne.*vex.*vez);
%     pxxe=mass(ele)*(pxxe-vex.*jex-jex.*vex+dne.*vex.*vex);
%     pxye=mass(ele)*(pxye-vex.*jey-jex.*vey+dne.*vex.*vey);
%     pyye=mass(ele)*(pyye-vey.*jey-jey.*vey+dne.*vey.*vey);
%     pyze=mass(ele)*(pyze-vey.*jez-jey.*vez+dne.*vey.*vez);
%     pzze=mass(ele)*(pzze-vez.*jez-jez.*vez+dne.*vez.*vez);

% I have no idea 

jix=jix*wpewce(1)*sqrt(mass(1));
jiy=jiy*wpewce(1)*sqrt(mass(1));
jiz=jiz*wpewce(1)*sqrt(mass(1));

vix = jix./dni;
viy = jiy./dni;
viz = jiz./dni;
vix(dni < 0.005) = 0;
viy(dni < 0.005) = 0;
viz(dni < 0.005) = 0;


jex=jex*wpewce(1)*sqrt(mass(1));
jey=jey*wpewce(1)*sqrt(mass(1));
jez=jez*wpewce(1)*sqrt(mass(1));

vex = jex./dne;
vey = jey./dne;
vez = jez./dne;

vex(dne < 0.005) = 0;
vey(dne < 0.005) = 0;
vez(dne < 0.005) = 0;


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

%     pi=(pxxi+pyyi+pzzi)/3;
%     pe=(pxxe+pyye+pzze)/3;
%     
%     
%     tmp=(jix.^2+jiy.^2+jiz.^2)/2;
%     eki = tmp./dni;
%     eki(dni < threshold) = 0;
% 
% 
%     tmp=(jex.^2+jey.^2+jez.^2)/2;
% %     dendiv,tmp,dne,eke,threshold
%     eke = tmp./dne;
%     eke(dni < threshold) = 0;
% 
%     tmp=(pxxi+pyyi+pzzi)/3;
% %     dendiv,tmp,dni,ti,threshold
%     ti = tmp./dni;
%     ti(dni < threshold) = 0;
% 
%     tmp=(pxxe+pyye+pzze)/3;
% %     dendiv,tmp,dne,te,threshold
%     te = tmp./dne;
%     te(dni < threshold) = 0;
% 
%      Vixz_o_divergence = 1;

%% Species 3,4
ion = 3;
ele = 4;

dni_h=squeeze(dns(:,:,ion(1)))*dfac(ion(1));
dne_h=squeeze(dns(:,:,ele(1)))*dfac(ele(1));

jix_h=squeeze(vxs(:,:,ion(1)))*dfac(ion(1));
jiy_h=squeeze(vys(:,:,ion(1)))*dfac(ion(1));
jiz_h=squeeze(vzs(:,:,ion(1)))*dfac(ion(1));

vix_h = jix_h./dni_h;
viy_h = jiy_h./dni_h;
viz_h = jiz_h./dni_h;
vix_h(dni_h < 0.005) = 0;
viy_h(dni_h < 0.005) = 0;
viz_h(dni_h < 0.005) = 0;

jex_h=squeeze(vxs(:,:,ele(1)))*dfac(ele(1)) ;
jey_h=squeeze(vys(:,:,ele(1)))*dfac(ele(1)) ;
jez_h=squeeze(vzs(:,:,ele(1)))*dfac(ele(1));

vex_h = jex_h./dne_h;
vey_h = jey_h./dne_h;
vez_h = jez_h./dne_h;

vex_h(dne_h < 0.005) = 0;
vey_h(dne_h < 0.005) = 0;
vez_h(dne_h < 0.005) = 0;

pxxi_h=squeeze(pxx(:,:,ion))*dfac(ion);
pyyi_h=squeeze(pyy(:,:,ion))*dfac(ion);
pzzi_h=squeeze(pzz(:,:,ion))*dfac(ion);
pxyi_h=squeeze(pxy(:,:,ion))*dfac(ion);
pxzi_h=squeeze(pxz(:,:,ion))*dfac(ion);
pyzi_h=squeeze(pyz(:,:,ion))*dfac(ion);

pxzi_h=mass(ion)*(pxzi_h-vix_h.*jiz_h);
pxxi_h=mass(ion)*(pxxi_h-vix_h.*jix_h);
pxyi_h=mass(ion)*(pxyi_h-vix_h.*jiy_h);
pyyi_h=mass(ion)*(pyyi_h-viy_h.*jiy_h);
pyzi_h=mass(ion)*(pyzi_h-viy_h.*jiz_h);
pzzi_h=mass(ion)*(pzzi_h-viz_h.*jiz_h);

pxxe_h=squeeze(pxx(:,:,ele))*dfac(ele);
pyye_h=squeeze(pyy(:,:,ele))*dfac(ele);
pzze_h=squeeze(pzz(:,:,ele))*dfac(ele);
pxye_h=squeeze(pxy(:,:,ele))*dfac(ele);
pxze_h=squeeze(pxz(:,:,ele))*dfac(ele);
pyze_h=squeeze(pyz(:,:,ele))*dfac(ele);

pxze_h=mass(ele)*(pxze_h-vex_h.*jez_h);
pxxe_h=mass(ele)*(pxxe_h-vex_h.*jex_h);
pxye_h=mass(ele)*(pxye_h-vex_h.*jey_h);
pyye_h=mass(ele)*(pyye_h-vey_h.*jey_h);
pyze_h=mass(ele)*(pyze_h-vey_h.*jez_h);
pzze_h=mass(ele)*(pzze_h-vez_h.*jez_h);



 jex_h=jex_h*wpewce(1)*sqrt(mass(1));
jey_h=jey_h*wpewce(1)*sqrt(mass(1));
jez_h=jez_h*wpewce(1)*sqrt(mass(1));

    vex_h = jex_h./dne_h;
    vey_h = jey_h./dne_h;
    vez_h = jez_h./dne_h;
    vex_h(dne_h == 0) = 0;
    vey_h(dne_h == 0) = 0;
    vez_h(dne_h == 0) = 0;

jix_h=jix_h*wpewce(1)*sqrt(mass(1));
jiy_h=jiy_h*wpewce(1)*sqrt(mass(1));
jiz_h=jiz_h*wpewce(1)*sqrt(mass(1));

    vix_h = jix_h./dni_h;
    viy_h = jiy_h./dni_h;
    viz_h = jiz_h./dni_h;
    vix_h(dni_h == 0) = 0;
    viy_h(dni_h == 0) = 0;
    viz_h(dni_h == 0) = 0;

pxxi_h=pxxi_h*wpewce(1)^2;
pyyi_h=pyyi_h*wpewce(1)^2;
pzzi_h=pzzi_h*wpewce(1)^2;
pxyi_h=pxyi_h*wpewce(1)^2;
pxzi_h=pxzi_h*wpewce(1)^2;
pyzi_h=pyzi_h*wpewce(1)^2;

pxxe_h=pxxe_h*wpewce(1)^2;
pyye_h=pyye_h*wpewce(1)^2;
pzze_h=pzze_h*wpewce(1)^2;
pxye_h=pxye_h*wpewce(1)^2;
pxze_h=pxze_h*wpewce(1)^2;
pyze_h=pyze_h*wpewce(1)^2;

%     pi_h=(pxxi_h+pyyi_h+pzzi_h)/3;
%     pe_h=(pxxe_h+pyye_h+pzze_h)/3;
%     
%     tmp=(jix_h.^2+jiy_h.^2+jiz_h.^2)/2;
%     eki_h = tmp./dni_h;
% %     dendiv,tmp,dni,eki,threshold
%     eki_h(dni_h < threshold) = 0;
% 
% 
%     tmp=(jex_h.^2+jey_h.^2+jez_h.^2)/2;
% %     dendiv,tmp,dne,eke,threshold
%     eke_h = tmp./dne_h;
%     eke_h(dni_h < threshold) = 0;
% 
%     tmp=(pxxi_h+pyyi_h+pzzi_h)/3;
% %     dendiv,tmp,dni,ti,threshold
%     ti_h = tmp./dni_h;
%     ti_h(dni_h < threshold) = 0;
% 
%     tmp=(pxxe_h+pyye_h+pzze_h)/3;
% %     dendiv,tmp,dne,te,threshold
%     te_h = tmp./dne_h;
%     te_h(dni_h < threshold) = 0;

%     Vixz_o_divergence

%     
%% Collect data in structures, r, b, e, ni, ne, vi, ve, ji, je, pi, pe
r.units = 'r/di';
r.x = xe;
r.z = ze;

b.units = 'B/B0';
b.x = bx;
b.y = by;
b.z = bz;
b.abs = sqrt(bx.^2 + by.^2 + bz.^2);

e.units = 'E/(vA*B0)';
e.x = ex;
e.y = ey;
e.z = ez;
e.abs = sqrt(ex.^2 + ey.^2 + ez.^2);
e.par = (ex.*bx + ey.*by + ez.*bz)./b.abs;
e.perp.x = ex - e.par.*bx./b.abs;
e.perp.y = ey - e.par.*by./b.abs;
e.perp.z = ez - e.par.*bz./b.abs;

% densoty
ni1 = dni;
ni2 = dni_h;

ne1 = dne;
ne2 = dne_h;

% velocity
ve1.units = 've/vA';
ve1.x = vex;
ve1.y = vey;
ve1.z = vez;
ve1.par = (ve1.x.*bx + ve1.y.*by + ve1.z.*bz)./b.abs;

ve2.units = 've/vA';
ve2.x = vex_h;
ve2.y = vey_h;
ve2.z = vez_h;
ve2.par = (ve2.x.*bx + ve2.y.*by + ve2.z.*bz)./b.abs;

vi1.units = 'vi/vA';
vi1.x = vix;
vi1.y = viy;
vi1.z = viz;
vi1.par = (vi1.x.*bx + vi1.y.*by + vi1.z.*bz)./b.abs;

vi2.units = 'vi/vA';
vi2.x = vix_h;
vi2.y = viy_h;
vi2.z = viz_h;
vi2.par = (vi2.x.*bx + vi2.y.*by + vi2.z.*bz)./b.abs;

% current
je1.units = 'je/(n0*vA)';
je1.x = jex;
je1.y = jex;
je1.z = jey;
je1.par = (je1.x.*bx + je1.y.*by + je1.z.*bz)./b.abs;

je2.units = 'je/(n0*vA)';
je2.x = jex_h;
je2.y = jex_h;
je2.z = jey_h;
je2.par = (je2.x.*bx + je2.y.*by + je2.z.*bz)./b.abs;

ji1.units = 'ji/(n0*vA)';
ji1.x = jix;
ji1.y = jiy;
ji1.z = jiz;
ji1.par = (ji1.x.*bx + ji1.y.*by + ji1.z.*bz)./b.abs;

ji2.units = 'ji/(n0*vA)';
ji2.x = jix_h;
ji2.y = jiy_h;
ji2.z = jiz_h;
ji2.par = (ji2.x.*bx + ji2.y.*by + ji2.z.*bz)./b.abs;

% pressure
pe1.units = 'pe/(...)';
pe1.xx = pxxe;
pe1.yy = pyye;
pe1.zz = pzze;
pe1.xy = pxye;
pe1.xz = pxze;
pe1.yz = pyze;
pe1.scalar = (pxxe + pyye + pzze)/3;

pe2.units = 'pe/(...)';
pe2.xx = pxxe_h;
pe2.yy = pyye_h;
pe2.zz = pzze_h;
pe2.xy = pxye_h;
pe2.xz = pxze_h;
pe2.yz = pyze_h;
pe2.scalar = (pxxe_h + pyye_h + pzze_h)/3;


pi1.units = 'pi/(...)';
pi1.xx = pxxi;
pi1.yy = pyyi;
pi1.zz = pzzi;
pi1.xy = pxyi;
pi1.xz = pxzi;
pi1.yz = pyzi;
pi1.scalar = (pxxi + pyyi + pzzi)/3;

pi2.units = 'pi/(...)';
pi2.xx = pxxi_h;
pi2.yy = pyyi_h;
pi2.zz = pzzi_h;
pi2.xy = pxyi_h;
pi2.xz = pxzi_h;
pi2.yz = pyzi_h;
pi2.scalar = (pxxi_h + pyyi_h + pzzi_h)/3;

% temperature
te1.units = 'te/(...)';
te1.xx = pe1.xx./ne1;
te1.yy = pe1.yy./ne1;
te1.zz = pe1.zz./ne1;
te1.xy = pe1.xy./ne1;
te1.xz = pe1.xz./ne1;
te1.yz = pe1.yz./ne1;
te1.scalar = (te1.xx + te1.yy + te1.zz)/3;

te2.units = 'te/(...)';
te2.xx = pe2.xx./ne2;
te2.yy = pe2.yy./ne2;
te2.zz = pe2.zz./ne2;
te2.xy = pe2.xy./ne2;
te2.xz = pe2.xz./ne2;
te2.yz = pe2.yz./ne2;
te2.scalar = (te2.xx + te2.yy + te2.zz)/3;

ti1.units = 'ti/(...)';
ti1.xx = pi1.xx./ni1;
ti1.yy = pi1.yy./ni1;
ti1.zz = pi1.zz./ni1;
ti1.xy = pi1.xy./ni1;
ti1.xz = pi1.xz./ni1;
ti1.yz = pi1.yz./ni1;
ti1.scalar = (ti1.xx + ti1.yy + ti1.zz)/3;

ti2.units = 'ti/(...)';
ti2.xx = pi2.xx./ni2;
ti2.yy = pi2.yy./ni2;
ti2.zz = pi2.zz./ni2;
ti2.xy = pi2.xy./ni2;
ti2.xz = pi2.xz./ni2;
ti2.yz = pi2.yz./ni2;
ti2.scalar = (ti2.xx + ti2.yy + ti2.zz)/3;
