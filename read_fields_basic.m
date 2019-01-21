function [time,r,e,b,ni,ne,ve,vi,je,ji,pe,pi,...
    dfac,teti,nnx,nnz, ...    
    wpewce,mass,it,dt,xmax,zmax,q] ...
    = read_fields_basic(txtfile,varargin)  
  % READ_FILES
  %   Reads files and return normalized quantities
  % [xe,ze,ex,ey,ez,bx,by,bz,ni,ne,ve,vi,je,ji,pe,pi,...
  %   dfac,teti,nnx,nnz, ...    
  %   wpewce,mass,it,dt,xmax,zmax,q] ...
  %   = READ_FILES_BASIC(pathFile,varargin)
  
  % defaults
  nss = 4; % numer of species
  numberel = 101; % number of bins for particle distributions
  groups = {[1 3],[2 4]}; %, i, e
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
      case 'group' % number of bins for particle distributions
        l = 2;
        groups = args{2};  
    end
    args = args((l+1):end);
    if isempty(args), break, end
  end
  
  nGroups = numel(groups);
  
  %% load data, data normalized units
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
  wpewce = fread(fid,1,'real*4');                               % wpe/wce 
  dfac = fread(fid,nss,'real*4');                               % dfac
  pxx_tmp = zeros(nnx,nnz,nss);
  pyy_tmp = zeros(nnx,nnz,nss);
  pzz_tmp = zeros(nnx,nnz,nss);
  pxy_tmp = zeros(nnx,nnz,nss);
  pxz_tmp = zeros(nnx,nnz,nss);
  pyz_tmp = zeros(nnx,nnz,nss);

  for is = 1:nss, pxx_tmp(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxx 
  for is = 1:nss, pyy_tmp(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyy 
  for is = 1:nss, pzz_tmp(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pzz 
  for is = 1:nss, pxy_tmp(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxy 
  for is = 1:nss, pxz_tmp(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pxz 
  for is = 1:nss, pyz_tmp(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % pyz
  remainder = fread(fid); 
 
  %% normalize 

  % fields
  c = 1;
  cwpe = sqrt(mass(2));
  cwpi = sqrt(mass(1));
  vAe = wpewce/c;
  vAp = vAe*sqrt(mass(1))/sqrt(mass(2));
  memi = mass(2)/mass(1);
  
  time = time*memi*wpewce^-1;
  xe=xe*cwpe/cwpi; % x/de -> x/di
  ze=ze*cwpe/cwpi; % z/de -> z/di
  r.units = 'r/di';
  r.x = xe;
  r.z = ze;
  bx=bx*wpewce; % B/(wpewce*B0) -> B/B0
  by=by*wpewce;
  bz=bz*wpewce;
  b.units = 'B/B0';
  b.x = bx;
  b.y = by;
  b.z = bz;
  ex=ex*sqrt(mass(1)/mass(2))*wpewce^2;
  ey=ey*sqrt(mass(1)/mass(2))*wpewce^2;
  ez=ez*sqrt(mass(1)/mass(2))*wpewce^2;
  e.units = 'EvAB0';
  e.x = ex;
  e.y = ey;
  e.z = ez;
  % moments
  for iSpecies = 1:nss
    dn(:,:,iSpecies) = dns(:,:,iSpecies)*dfac(iSpecies);
    jx(:,:,iSpecies) = q(iSpecies)*vxs(:,:,iSpecies)*dfac(iSpecies)*wpewce*sqrt(1/memi);
    jy(:,:,iSpecies) = q(iSpecies)*vys(:,:,iSpecies)*dfac(iSpecies)*wpewce*sqrt(1/memi);
    jz(:,:,iSpecies) = q(iSpecies)*vzs(:,:,iSpecies)*dfac(iSpecies)*wpewce*sqrt(1/memi);
    
    vx(:,:,iSpecies) = jx(:,:,iSpecies)./dn(:,:,iSpecies)/q(iSpecies);
    vy(:,:,iSpecies) = jy(:,:,iSpecies)./dn(:,:,iSpecies)/q(iSpecies);
    vz(:,:,iSpecies) = jz(:,:,iSpecies)./dn(:,:,iSpecies)/q(iSpecies);
    vx(dn==0) = 0;
    vy(dn==0) = 0;
    vz(dn==0) = 0;
    
    if 1
    pxx(:,:,iSpecies) = memi^-1*wpewce^2*( pxx_tmp(:,:,iSpecies)*dfac(iSpecies) - vx(:,:,iSpecies).*jx(:,:,iSpecies)/q(iSpecies) );    
    pyy(:,:,iSpecies) = memi^-1*wpewce^2*( pyy_tmp(:,:,iSpecies)*dfac(iSpecies) - vy(:,:,iSpecies).*jy(:,:,iSpecies)/q(iSpecies) );
    pzz(:,:,iSpecies) = memi^-1*wpewce^2*( pzz_tmp(:,:,iSpecies)*dfac(iSpecies) - vz(:,:,iSpecies).*jz(:,:,iSpecies)/q(iSpecies) );
    pxy(:,:,iSpecies) = memi^-1*wpewce^2*( pxy_tmp(:,:,iSpecies)*dfac(iSpecies) - vx(:,:,iSpecies).*jy(:,:,iSpecies)/q(iSpecies) );
    pxz(:,:,iSpecies) = memi^-1*wpewce^2*( pxz_tmp(:,:,iSpecies)*dfac(iSpecies) - vx(:,:,iSpecies).*jz(:,:,iSpecies)/q(iSpecies) );
    pyz(:,:,iSpecies) = memi^-1*wpewce^2*( pyz_tmp(:,:,iSpecies)*dfac(iSpecies) - vy(:,:,iSpecies).*jz(:,:,iSpecies)/q(iSpecies) );
    else
    pxx(:,:,iSpecies) = memi^-1*wpewce^2*( pxx_tmp(:,:,iSpecies) - vxs(:,:,iSpecies).*vxs(:,:,iSpecies) )*dfac(iSpecies)^1;    
    pyy(:,:,iSpecies) = memi^-1*wpewce^2*( pyy_tmp(:,:,iSpecies) - vys(:,:,iSpecies).*vys(:,:,iSpecies) )*dfac(iSpecies)^1;
    pzz(:,:,iSpecies) = memi^-1*wpewce^2*( pzz_tmp(:,:,iSpecies) - vzs(:,:,iSpecies).*vzs(:,:,iSpecies) )*dfac(iSpecies)^1;
    pxy(:,:,iSpecies) = memi^-1*wpewce^2*( pxy_tmp(:,:,iSpecies) - vxs(:,:,iSpecies).*vys(:,:,iSpecies) )*dfac(iSpecies)^1;
    pxz(:,:,iSpecies) = memi^-1*wpewce^2*( pxz_tmp(:,:,iSpecies) - vxs(:,:,iSpecies).*vzs(:,:,iSpecies) )*dfac(iSpecies)^1;
    pyz(:,:,iSpecies) = memi^-1*wpewce^2*( pyz_tmp(:,:,iSpecies) - vys(:,:,iSpecies).*vzs(:,:,iSpecies) )*dfac(iSpecies)^1;
    end
    
    if 0
      %%
      A = pzz_tmp(:,:,iSpecies)*dfac(iSpecies);
      B1 = vz(:,:,iSpecies);
      B2 = jz(:,:,iSpecies);
      B = B1.*B2;
      C = mass(iSpecies)*wpewce^2;
      P = C*(A-B);
      nrows = 3;
      ncols = 2;
      npanels = nrows*ncols;
      for ipanel = 1:npanels
        h(ipanel) = subplot(nrows,ncols,ipanel);
      end

      isub = 1;
      hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,A'); colorbar('peer',hca);
      hca.Title.String = 'pzz tmp * dfac';
      
      hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,B1'); colorbar('peer',hca);
      hca.Title.String = 'vz';
      
      hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,B2'); colorbar('peer',hca);
      hca.Title.String = 'jz';
      
      hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,B'); colorbar('peer',hca);
      hca.Title.String = 'vz*jz';   
      
      hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,C'); colorbar('peer',hca);
      hca.Title.String = 'm*(wpe/wce)^2';
      
      hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,P'); colorbar('peer',hca);
      hca.Title.String = 'm*(wpe/wce)^2 * (pzz tmp * dfac - vz*jz)';
      
    end
  end
  if 0
      %%   
      nrows = 4;
      ncols = 5;
      npanels = nrows*ncols;
      for ipanel = 1:npanels
        h(ipanel) = subplot(nrows,ncols,ipanel);
      end

      isub = 1;
      for iSpecies = 1:4
        hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,pzz_tmp(:,:,iSpecies)'); colorbar('peer',hca);
        hca.Title.String = sprintf('pzz code, species %d',iSpecies);
        colorbar('peer',hca);  
        hca.CLim = max(abs(hca.CLim))*[-1 1]; 
        
        hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,vzs(:,:,iSpecies)'.*vzs(:,:,iSpecies)'); colorbar('peer',hca);
        hca.Title.String = sprintf('vz*vz code, species %d',iSpecies);
        colorbar('peer',hca);  
        hca.CLim = max(abs(hca.CLim))*[-1 1]; 
        
        hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,dfac(iSpecies)*pzz_tmp(:,:,iSpecies)'); colorbar('peer',hca);
        hca.Title.String = sprintf('dfac * pzz code, species %d',iSpecies);
        colorbar('peer',hca);  
        hca.CLim = max(abs(hca.CLim))*[-1 1]; 
        
        hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,vz(:,:,iSpecies).*jz(:,:,iSpecies)'); colorbar('peer',hca);
        hca.Title.String = sprintf('vz*jz, species %d',iSpecies);
        colorbar('peer',hca);  
        hca.CLim = max(abs(hca.CLim))*[-1 1]; 
        hca.CLim = 0.1*[-1 1];
        
        hca = h(isub); isub = isub + 1; imagesc(hca,xe,ze,pzz(:,:,iSpecies)'); colorbar('peer',hca);
        hca.Title.String = sprintf('pzz, species %d',iSpecies);
        colorbar('peer',hca);  
        hca.CLim = max(abs(hca.CLim))*[-1 1];
        
        colormap(cn.cmap('blue_red'))
      end
      
      
    end
  %vx = jx./dn;
  %vy = jy./dn;
  %vz = jz./dn;
  %vx(dn == 0) = 0;
  %vy(dn == 0) = 0;
  %vz(dn == 0) = 0;   
  
  %% sum over species, and collect for output
  
    for iGroup = 1:nGroups
      membersGroup = groups{iGroup};
      nMembers = numel(membersGroup); 
        if numel(unique(q(membersGroup))) == 1
          q_species = q(membersGroup(1));
        else
          error('Something wrong with species charge.');
        end
      if mass(groups{iGroup}) == max(unique(mass)) % ions   
        % total
        % n      
        ni.tot = squeeze(sum(dn(:,:,membersGroup),3));
        % j
        ji.x = squeeze(sum(jx(:,:,membersGroup),3));
        ji.y = squeeze(sum(jy(:,:,membersGroup),3));
        ji.z = squeeze(sum(jz(:,:,membersGroup),3));
        % v
        if 0 % can not be summed like this, because dfac is different
          vi.x = squeeze(sum(vx(:,:,membersGroup),3));
          vi.y = squeeze(sum(vy(:,:,membersGroup),3));
          vi.z = squeeze(sum(vz(:,:,membersGroup),3));
        else % get total velocity from total current: j/n
          vi.x = ji.x./ni.tot/q_species;
          vi.y = ji.y./ni.tot/q_species;
          vi.z = ji.z./ni.tot/q_species;
        end
        % p
        pi.xx = squeeze(sum(pxx(:,:,membersGroup),3));
        pi.yy = squeeze(sum(pyy(:,:,membersGroup),3));
        pi.zz = squeeze(sum(pzz(:,:,membersGroup),3));
        pi.xy = squeeze(sum(pxy(:,:,membersGroup),3));
        pi.xz = squeeze(sum(pxz(:,:,membersGroup),3));
        pi.yz = squeeze(sum(pyz(:,:,membersGroup),3));
        pi.scalar = (pi.xx + pi.yy + pi.zz)/3;
        % for each population separately
        for iMember = 1:nMembers
          c_eval('ni.s! = squeeze(dn(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('vi.x! = squeeze(vx(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('vi.y! = squeeze(vy(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('vi.x! = squeeze(vz(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('ji.x! = squeeze(jx(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('ji.y! = squeeze(jy(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('ji.z! = squeeze(jz(:,:,groups{?}(!)));',iGroup,iMember)  
          c_eval('pi.xx! = squeeze(pxx(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pi.yy! = squeeze(pyy(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pi.zz! = squeeze(pzz(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pi.xy! = squeeze(pxy(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pi.xz! = squeeze(pxz(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pi.yz! = squeeze(pyz(:,:,groups{?}(!)));',iGroup,iMember)       
          c_eval('pi.scalar? = (pi.xx? + pi.yy? + pi.zz?)/3;',iMember)
        end
      elseif mass(groups{iGroup}) == min(unique(mass)) % electrons   
        % total
        % n      
        ne.tot = squeeze(sum(dn(:,:,membersGroup),3));
        ve.x = squeeze(sum(vx(:,:,membersGroup),3));
        ve.y = squeeze(sum(vy(:,:,membersGroup),3));
        ve.z = squeeze(sum(vz(:,:,membersGroup),3));
        % j
        je.x = squeeze(sum(jx(:,:,membersGroup),3));
        je.y = squeeze(sum(jy(:,:,membersGroup),3));
        je.z = squeeze(sum(jz(:,:,membersGroup),3));        
        % v
        if 0 % can not be summed like this, because dfac is different
          ve.x = squeeze(sum(vx(:,:,membersGroup),3));
          ve.y = squeeze(sum(vy(:,:,membersGroup),3));
          ve.z = squeeze(sum(vz(:,:,membersGroup),3));
        else % get total velocity from total current: j/n
          ve.x = je.x./ni.tot/q_species;
          ve.y = je.y./ni.tot/q_species;
          ve.z = je.z./ni.tot/q_species;
        end
        
        % p
        pe.xx = squeeze(sum(pxx(:,:,membersGroup),3));
        pe.yy = squeeze(sum(pyy(:,:,membersGroup),3));
        pe.zz = squeeze(sum(pzz(:,:,membersGroup),3));
        pe.xy = squeeze(sum(pxy(:,:,membersGroup),3));
        pe.xz = squeeze(sum(pxz(:,:,membersGroup),3));
        pe.yz = squeeze(sum(pyz(:,:,membersGroup),3));
        pe.scalar = (pe.xx + pe.yy + pe.zz)/3;
        % for each population separately
        for iMember = 1:nMembers
          c_eval('ne.s! = squeeze(dn(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('ve.x! = squeeze(vx(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('ve.y! = squeeze(vy(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('ve.x! = squeeze(vz(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('je.x! = squeeze(jx(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('je.y! = squeeze(jy(:,:,groups{?}(!)));',iGroup,iMember)
          c_eval('je.z! = squeeze(jz(:,:,groups{?}(!)));',iGroup,iMember)  
          c_eval('pe.xx! = squeeze(pxx(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pe.yy! = squeeze(pyy(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pe.zz! = squeeze(pzz(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pe.xy! = squeeze(pxy(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pe.xz! = squeeze(pxz(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pe.yz! = squeeze(pyz(:,:,groups{?}(!)));',iGroup,iMember) 
          c_eval('pe.scalar? = (pe.xx? + pe.yy? + pe.zz?)/3;',iMember)
        end
      end    
    end
