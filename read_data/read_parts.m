function [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
    = read_parts(txtfile,varargin)
  % READ_PARTS
  % [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = READ_PARTS(txtfile,nss)
  % nss = number of species, if not given, use default 4
  
  % defaults
  nss = 4; % numer of species
  numberel = 101; % number of bins for particle distributions
  
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
    end
    args = args((l+1):end);
    if isempty(args), break, end
  end
  
  fid = fopen(txtfile,'r','ieee-le');

  % preallocate
  fxy = zeros(numberel,numberel,nss);
  fxz = zeros(numberel,numberel,nss);
  fyz = zeros(numberel,numberel,nss);
  axes = zeros(numberel,nss);

  % read file
  fread(fid,1,'integer*8');
  for is = 1:nss, axes(:,is) = fread(fid,numberel,'real*4'); end
  xlo = fread(fid,1,'real*4');                                 
  xhi = fread(fid,1,'real*4');                                 
  zlo = fread(fid,1,'real*4');
  zhi = fread(fid,1,'real*4');
  ic = fread(fid,nss,'integer*4'); % long
  fxyz = fread(fid,numberel*numberel*numberel*nss,'real*4');
  fxyz = reshape(fxyz,numberel,numberel,numberel,nss);
  for is = 1:nss, fxy(:,:,is) = fread(fid,[numberel numberel],'real*4'); end  % vxs 
  for is = 1:nss, fxz(:,:,is) = fread(fid,[numberel numberel],'real*4'); end  % vxs 
  for is = 1:nss, fyz(:,:,is) = fread(fid,[numberel numberel],'real*4'); end  % vxs 
  vxa = fread(fid,nss,'real*4');
  vya = fread(fid,nss,'real*4');
  vza = fread(fid,nss,'real*4');

  fclose(fid);