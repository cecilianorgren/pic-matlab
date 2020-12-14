function h5write_fields_complement(pic_orig)
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

timesteps = pic_orig.twpe;
for it = 1:pic_orig.nt
  pic_tmp = pic_orig(it);
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB = sum(Babs(:).^2)/2;           
  A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);        
  [Ainds,Avals] = saddle(A,'sort');
  xXline = pic_tmp.xi(Ainds(1,1));
  zXline = pic_tmp.zi(Ainds(1,2));  
  AXline = Avals(1);
  EyXline = mean(mean(pic_tmp.xlim(xXline+[-0.1 0.1]).zlim(zXline+[-0.1 0.1]).Ey));
  
  % Write data to file
  h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)
  h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'Axline',AXline)
  h5write_attr(pic_tmp,pic_tmp.twci,'xline_position',[xXline' zXline'])
  h5write_attr(pic_tmp,pic_tmp.twci,'RE',EyXline)  
  h5write_attr(pic_tmp,pic_tmp.twci,'UB',UB)
  
end
