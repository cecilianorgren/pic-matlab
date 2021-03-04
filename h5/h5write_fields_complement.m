function h5write_fields_complement(pic_orig,varargin)
% H5WRITE_FIELDS_COMPLEMENT Writes complementory info to h5 fields file. 
%
%   H5WRITE_FIELDS_COMPLEMENT(picobject)
%
%   Data added:
%     A - magnetic vector potential Ay
%     UB - Total magnetic energy
%     A_xline - A value at main X line
%     x_xline - x location of main X line
%     z_xline - z location of main X line
%     Ey_xline - Ey at main X line (reconnection electric field)
%

doGroup = 0;
if nargin > 1 && strcmp(varargin{1},'group') % additional input given
  doGroup = 1;
  groupSpecies = varargin{2};
  nGroups = numel(groupSpecies);
end

doGroup = 0;
groups = {[1],[2],[3 5],[4 6]};
groupNames = {'i_hot','e_hot','i_cold','e_cold'};
  
timesteps = pic_orig.twpe;
for it = 1:pic_orig.nt
  pic_tmp = pic_orig(it);
  if 1
    Bx = pic_tmp.Bx;
    By = pic_tmp.By;
    Bz = pic_tmp.Bz;
    Babs = sqrt(Bx.^2+By.^2+Bz.^2);
    % Magnetic field energy density
    UB = sum(Babs(:).^2)/2;
    
    % Write data to file
    h5write_attr(pic_tmp,pic_tmp.twci,'UB',UB)
  end
  % Plasma energy densities
  if 1
    ne = pic_tmp.ne;  
    pe = pic_tmp.pe;
    vex = pic_tmp.vex;
    vey = pic_tmp.vey;
    vez = pic_tmp.vez;
    ni = pic_tmp.ni;
    pi = pic_tmp.pi;
    vix = pic_tmp.vix;
    viy = pic_tmp.viy;
    viz = pic_tmp.viz;
    Uke = pic_tmp.mass(2)/pic_tmp.mass(1)*0.5*ne.*(vex.^2 + vey.^2 + vez.^2);
    Uki = pic_tmp.mass(1)/pic_tmp.mass(1)*0.5*ni.*(vix.^2 + viy.^2 + viz.^2);
    Ute = 3/2*pe;
    Uti = 3/2*pi;

    Uke = nansum(Uke(:));
    Ute = nansum(Ute(:));
    Uki = nansum(Uki(:));
    Uti = nansum(Uti(:));
    
    % Write data to file
    h5write_attr(pic_tmp,pic_tmp.twci,'Ute',Ute)
    h5write_attr(pic_tmp,pic_tmp.twci,'Uti',Uti)
    h5write_attr(pic_tmp,pic_tmp.twci,'Uke',Uke)
    h5write_attr(pic_tmp,pic_tmp.twci,'Uki',Uki)
  end
  
  if doGroup % NOT IMPLEMENTED
    for iGroup = 1:nGroups
      iSpecies = groupSpecies{iGroup};
      n = pic_tmp.n(iSpecies);
      p = pic_tmp.p(iSpecies);
      vx = pic_tmp.vx(iSpecies);
      vy = pic_tmp.vy(iSpecies);
      vz = pic_tmp.vz(iSpecies);      
      Uk = pic_tmp.mass(iSpecies(1))/pic_tmp.mass(1)*0.5*n.*(vx.^2 + vy.^2 + vz.^2);
      Ut = 3/2*p.scalar;         
      h5write_attr(pic_tmp,pic_tmp.twci,sprintf('Ut%s',groupNames{iGroup}),Ut)
      h5write_attr(pic_tmp,pic_tmp.twci,sprintf('Uk%s',groupNames{iGroup}),Uk)      
    end
  end
  
  % X line and reconnection rate
  if 1
    A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);        
    [Ainds,Avals] = saddle(A,'sort');
    xXline = pic_tmp.xi(Ainds(1,1));
    zXline = pic_tmp.zi(Ainds(1,2));  
    AXline = Avals(1);
    EyXline = mean(mean(pic_tmp.xlim(xXline+[-0.1 0.1]).zlim(zXline+[-0.1 0.1]).Ey));
    % Write data to file
    h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)
    h5write_attr(pic_tmp,pic_tmp.twci,'Axline',AXline)
    h5write_attr(pic_tmp,pic_tmp.twci,'xline_position',[xXline' zXline'])
    h5write_attr(pic_tmp,pic_tmp.twci,'RE',EyXline)  
  end
  
end
