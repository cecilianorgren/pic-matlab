% check 


sim_tmp = df04.ilim(40);%.zlim([-4 4]);
species = [3 5];

% combined
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(species);    
pdyn = sim_tmp.mass(species(1))/sim_tmp.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
p = (pxx+pyy+pzz)/3; % scalar pressure
UT_comb = 3/2*nansum(p(:));
UK_comb = nansum(pdyn(:));

clear UT_sep UK_sep

for iSpecies = 1:numel(species)  
  [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(species(iSpecies));    
  pdyn = sim_tmp.mass(species(iSpecies))/sim_tmp.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
  p = (pxx+pyy+pzz)/3; % scalar pressure
  UT_sep(iSpecies) = 3/2*nansum(p(:));
  UK_sep(iSpecies) = nansum(pdyn(:));
end

[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv(species); 
vv = (vxx+vyy+vzz)/2;
vvsum = sum(vv(:)) % "right answer"

sum_sep = sum(UT_sep) + sum(UK_sep)
sum_comb = UT_comb + UK_comb
%%
[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv(1); 
vv = (vxx+vyy+vzz)/2;
% these are the same
vvsum = sum(vv(:))
UT_sep(1)+UK_sep(1)

%%
sim_tmp = df04.ilim(30).zlim([-1 1]);
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp([1 3]);    
pdyn = sim_tmp.mass(1)/sim_tmp.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
p = (pxx+pyy+pzz)/3; % scalar pressure
UT_comb = 3/2*nansum(p(:));
UK_comb = nansum(pdyn(:));

[vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([1 3]); % this addition is straightforward
vv = (vxx+vyy+vzz)/2; 

% these are the same now
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

%% Plot current timeseries

plot(df04.twci,df04.UK(135)+df04.UT(135),df04.twci,df04.UK(1)+df04.UT(1)+df04.UK(3)+df04.UT(3)+df04.UK(5)+df04.UT(5),'--')
legend('combined','separate and added')

%%
plot(df04.twci,df04.UK(35)+df04.UT(35),df04.twci,df04.UK(3)+df04.UT(3)+df04.UK(5)+df04.UT(5),'--')
legend('combined','separate and added')

%% Check pressure balance
ii = 10; % time indice
pic = df04.i(ii,100:200,:);
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = pic.njp(1);
vzz = pic.vzz(1);
vyy = pic.vyy(1);
vv_diag = pic.vv_diag(1);
PBx = pic.Bx.^2/2;
plot(pic.zi,squeeze(mean(pzz,2)),'-.',...
  pic.zi,squeeze(mean(PBx,2)),...
  pic.zi,squeeze(mean(vzz,2)),...
  pic.zi,squeeze(mean(vv_diag,2)),...
  pic.zi,squeeze(mean(vyy,2)))
legend({'pzz','B^2/2','vzz','(vxx+vyy+vzz)/3','vyy','vxx'})
%% Check vv
ii = 10; % time indice
pic = df04.i(ii,100:200,:);
[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = pic.njp(1);
vxx = pic.vxx(1);
vyy = pic.vyy(1);
vzz = pic.vzz(1);
vv_diag = pic.vv_diag(1);
PBx = pic.Bx.^2/2;
plot(...
  pic.zi,squeeze(mean(vxx,2)),...
  pic.zi,squeeze(mean(vyy,2)),...
  pic.zi,squeeze(mean(vzz,2)),...
  pic.zi,squeeze(mean(vv_diag,2)))
legend({'vxx','vyy','vzz','(vxx+vyy+vzz)/3'})
%%
plot(df04.twci,df04.UK(46)+df04.UT(46),df04.twci,df04.UK(4)+df04.UT(4)+df04.UK(6)+df04.UT(6),'--')
legend('combined','separate and added')

%%
plot(df08.twci,df08.UK(13)+df08.UT(13),df08.twci,df08.UK(1)+df08.UT(1)+df08.UK(3)+df08.UT(3),'--')
legend('combined','separate and added')

%%
plot(df08.twci,df08.UK(24)+df08.UT(24),df08.twci,df08.UK(2)+df08.UT(2)+df08.UK(4)+df08.UT(4),'--')
legend('combined','separate and added')