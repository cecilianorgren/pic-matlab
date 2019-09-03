function varargout = read_fields_adaptive(txtfile,varargin)  
%[xe,ze,e,b,ni1,ne1,ni2,ne2,vi1,ve1,vi2,ve2,ji1,je1,ji2,je2,...
    %pi1,pe1,pi2,pe2,ti1,te1,ti2,te2,dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] ...

    % READ_FILES
  %   Reads files and return normalized quantities
 
  %% defaults
  nss = 6; % number of species
  doFixNegDensities = 0;
  %groups = {[1 3],[2 4]}; % electron and ions
  
  % read input
  nargs = numel(varargin);
  args = varargin;
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
      case 'rem_neg_n' % number of bins for particle distributions
        l = 1;
        doFixNegDensities = 1;
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
  vxx = zeros(nnx,nnz,nss); % vxvx
  vyy = zeros(nnx,nnz,nss); % vyvy
  vzz = zeros(nnx,nnz,nss); % vzvz
  vxy = zeros(nnx,nnz,nss); % vxvy
  vxz = zeros(nnx,nnz,nss); % vxvz
  vyz = zeros(nnx,nnz,nss); % vyvz

  for is = 1:nss, vxx(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxx 
  for is = 1:nss, vyy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vyy 
  for is = 1:nss, vzz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vzz 
  for is = 1:nss, vxy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxy 
  for is = 1:nss, vxz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxz 
  for is = 1:nss, vyz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vyz
  remainder = fread(fid);
  st = fclose(fid);
  
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
  xe = xe/sqrt(mass(1));
  ze = ze/sqrt(mass(1));
  bx = bx*wpewce(1);
  by = by*wpewce(1);
  bz = bz*wpewce(1);
  ex = ex*sqrt(mass(1))*wpewce(1)^2;
  ey = ey*sqrt(mass(1))*wpewce(1)^2;
  ez = ez*sqrt(mass(1))*wpewce(1)^2;
  
  time = time/wpewce/mass(1);
  
  % moments
  % initialize matrices
  n = zeros(numel(xe),numel(ze),nss);
  jx = zeros(numel(xe),numel(ze),nss);
  jy = zeros(numel(xe),numel(ze),nss);
  jz = zeros(numel(xe),numel(ze),nss);
  vx = zeros(numel(xe),numel(ze),nss);
  vy = zeros(numel(xe),numel(ze),nss);
  vz = zeros(numel(xe),numel(ze),nss);
  pxx = zeros(numel(xe),numel(ze),nss);
  pyy = zeros(numel(xe),numel(ze),nss);
  pzz = zeros(numel(xe),numel(ze),nss);
  pxy = zeros(numel(xe),numel(ze),nss);
  pxz = zeros(numel(xe),numel(ze),nss);
  pyz = zeros(numel(xe),numel(ze),nss);
  txx = zeros(numel(xe),numel(ze),nss);
  tyy = zeros(numel(xe),numel(ze),nss);
  tzz = zeros(numel(xe),numel(ze),nss);
  txy = zeros(numel(xe),numel(ze),nss);
  txz = zeros(numel(xe),numel(ze),nss);
  tyz = zeros(numel(xe),numel(ze),nss);
  
  for iSpecies = 1:nss
    % density, n
    n(:,:,iSpecies) = dns(:,:,iSpecies)*dfac(iSpecies);
    % flux, j = nv
    jx(:,:,iSpecies) = vxs(:,:,iSpecies)*dfac(iSpecies)*wpewce(1)*sqrt(mass(1)); % c/vA0
    jy(:,:,iSpecies) = vys(:,:,iSpecies)*dfac(iSpecies)*wpewce(1)*sqrt(mass(1));
    jz(:,:,iSpecies) = vzs(:,:,iSpecies)*dfac(iSpecies)*wpewce(1)*sqrt(mass(1));
    % velocity, v
    vx(:,:,iSpecies) = jx(:,:,iSpecies)./n(:,:,iSpecies);
    vy(:,:,iSpecies) = jy(:,:,iSpecies)./n(:,:,iSpecies);
    vz(:,:,iSpecies) = jz(:,:,iSpecies)./n(:,:,iSpecies);    
    % do this here already to not mess up something with the pressure I think
    vx(n==0) = 0;
    vy(n==0) = 0;
    vz(n==0) = 0;
    % pressure
    pxx(:,:,iSpecies) = mass(iSpecies)*wpewce^2*( vxx(:,:,iSpecies)*dfac(iSpecies) - vxs(:,:,iSpecies).*vxs(:,:,iSpecies)./n(:,:,iSpecies)*dfac(iSpecies)^2 );
    pyy(:,:,iSpecies) = mass(iSpecies)*wpewce^2*( vyy(:,:,iSpecies)*dfac(iSpecies) - vys(:,:,iSpecies).*vys(:,:,iSpecies)./n(:,:,iSpecies)*dfac(iSpecies)^2 );
    pzz(:,:,iSpecies) = mass(iSpecies)*wpewce^2*( vzz(:,:,iSpecies)*dfac(iSpecies) - vzs(:,:,iSpecies).*vzs(:,:,iSpecies)./n(:,:,iSpecies)*dfac(iSpecies)^2 );
    pxy(:,:,iSpecies) = mass(iSpecies)*wpewce^2*( vxy(:,:,iSpecies)*dfac(iSpecies) - vxs(:,:,iSpecies).*vys(:,:,iSpecies)./n(:,:,iSpecies)*dfac(iSpecies)^2 );
    pxz(:,:,iSpecies) = mass(iSpecies)*wpewce^2*( vxz(:,:,iSpecies)*dfac(iSpecies) - vxs(:,:,iSpecies).*vzs(:,:,iSpecies)./n(:,:,iSpecies)*dfac(iSpecies)^2 );
    pyz(:,:,iSpecies) = mass(iSpecies)*wpewce^2*( vyz(:,:,iSpecies)*dfac(iSpecies) - vys(:,:,iSpecies).*vzs(:,:,iSpecies)./n(:,:,iSpecies)*dfac(iSpecies)^2 );
    % temperature
    txx(:,:,iSpecies) = pxx(:,:,iSpecies)./n(:,:,iSpecies);
    tyy(:,:,iSpecies) = pyy(:,:,iSpecies)./n(:,:,iSpecies);
    tzz(:,:,iSpecies) = pzz(:,:,iSpecies)./n(:,:,iSpecies);
    txy(:,:,iSpecies) = pxy(:,:,iSpecies)./n(:,:,iSpecies);
    txz(:,:,iSpecies) = pxz(:,:,iSpecies)./n(:,:,iSpecies);
    tyz(:,:,iSpecies) = pyz(:,:,iSpecies)./n(:,:,iSpecies);
  end
  for iSpecies = nss:-1:1
    mass(iSpecies) = mass(iSpecies)/mass(1);
  end
%% Fix data
if doFixNegDensities
  pxxe()
end

%% Collect data in structures, x, z, b, e, ni, ne, vi, ve, ji, je, pi, pe
x = xe;
z = ze;

B.x = bx;
B.y = by;
B.z = bz;

E.x = ex;
E.y = ey;
E.z = ez;

species_str = {};
for iPop = 1:nss/2
  species_str{end+1} = sprintf('i%g',iPop);
  species_str{end+1} = sprintf('e%g',iPop);
end
comp_tens = {'xx','xy','xz','yy','yz','zz'};
comp_vec = {'x','y','z'};
varstrs_species = {};
for iSpecies = 1:nss
  % scalar
  eval(sprintf('n%s = squeeze(n(:,:,%g));',species_str{iSpecies},iSpecies))
  % vectors
  for iComp = 1:numel(comp_vec)
    eval(sprintf('v%s.%s = squeeze(v%s(:,:,%g));',species_str{iSpecies},comp_vec{iComp},comp_vec{iComp},iSpecies))
    eval(sprintf('j%s.%s = squeeze(j%s(:,:,%g));',species_str{iSpecies},comp_vec{iComp},comp_vec{iComp},iSpecies))
  end
  % tensors
  for iComp = 1:numel(comp_tens)
    eval(sprintf('p%s.%s = squeeze(p%s(:,:,%g));',species_str{iSpecies},comp_tens{iComp},comp_tens{iComp},iSpecies))
    eval(sprintf('t%s.%s = squeeze(t%s(:,:,%g));',species_str{iSpecies},comp_tens{iComp},comp_tens{iComp},iSpecies))
  end
  % collect strings for varargout
  varstrs_species{end+1,1} = sprintf('n%s',species_str{iSpecies});
  varstrs_species{end+1,1} = sprintf('v%s',species_str{iSpecies});
  varstrs_species{end+1,1} = sprintf('j%s',species_str{iSpecies});
  varstrs_species{end+1,1} = sprintf('p%s',species_str{iSpecies});
  varstrs_species{end+1,1} = sprintf('t%s',species_str{iSpecies});
end

% Collect data in single varargout

varstrs_common = {'x','z','E','B','dfac','teti','nnx','nnz','wpewce','mass','it','time','dt','xmax','zmax','q'};

irow = 0;
for iout = 1:numel(varstrs_common)
  irow = irow + 1;
  varargout_{irow,1} = varstrs_common{iout};
  varargout_{irow,2} = eval(varstrs_common{iout});
end
for iout = 1:numel(varstrs_species)
  irow = irow + 1;
  varargout_{irow,1} = varstrs_species{iout};
  varargout_{irow,2} = eval(varstrs_species{iout});
end
varargout = {varargout_};
