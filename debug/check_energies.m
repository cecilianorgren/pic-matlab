% check 

% combined
sim_tmp = df04.ilim(30).zlim([-15 -10]);
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp([1 3 5]);    
pdyn = sim.mass(iSpecies)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
p = (pxx+pyy+pzz)/3; % scalar pressure
UT_comb = 3/2*nansum(p(:));
UK_comb = nansum(pdyn(:));

species = [1 3 5];
for iSpecies = 1:numel(species)  
  [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(species(iSpecies));    
  pdyn = sim.mass(iSpecies)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
  p = (pxx+pyy+pzz)/3; % scalar pressure
  UT_sep(iSpecies) = 3/2*nansum(p(:));
  UK_sep(iSpecies) = nansum(pdyn(:));
end

[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([1 3 5]); 
vv = (vxx+vyy+vzz)/2;
vvsum = sum(vv(:));

sum(UT_sep) + sum(UK_sep)
UT_comb + UK_comb
%%
[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv(1); 
vv = (vxx+vyy+vzz)/2;
% these are the same
vvsum = sum(vv(:))
UT_sep(1)+UK_sep(1)

%%
sim_tmp = df04.ilim(30).zlim([-1 1]);
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp([1 3 5]);    
pdyn = sim.mass(1)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
p = (pxx+pyy+pzz)/3; % scalar pressure
UT_comb = 3/2*nansum(p(:));
UK_comb = nansum(pdyn(:));

[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([1 3 5]); % this addition is straightforward
vv = (vxx+vyy+vzz)/2; 

% these are not the same
vvsum = sum(vv(:))
UT_comb+UK_comb
%%
sim_tmp = df04.ilim(30).zlim([-1 1]);
tic
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp([1 3 5]);
[n1,jx1,jy1,jz1,pxx1,pxy1,pxz1,pyy1,pyz1,pzz1] = sim_tmp.njp([1]);
[n3,jx3,jy3,jz3,pxx3,pxy3,pxz3,pyy3,pyz3,pzz3] = sim_tmp.njp([3]);
[n5,jx5,jy5,jz5,pxx5,pxy5,pxz5,pyy5,pyz5,pzz5] = sim_tmp.njp([5]);
[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([1 3 5]); 
[vxx1,vxy1,vxz1,vyy1,vyz1,vzz1] = sim_tmp.vv(1); 
[vxx3,vxy3,vxz3,vyy3,vyz3,vzz3] = sim_tmp.vv(3); 
[vxx5,vxy5,vxz5,vyy5,vyz5,vzz5] = sim_tmp.vv(5); 
toc

% these are the sme
sum(sum(n-n1-n3-n5))/sum(n(:))
sum(sum(jy-jy1-jy3-jy5))/sum(jy(:))
sum(sum(vyy-vyy1-vyy3-vyy5))/sum(vyy(:))