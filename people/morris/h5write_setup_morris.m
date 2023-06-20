if 0
xlim = mean(pic.xi)+0.2*[-1 1];
twpe = 0;
A = mean(pic.twpelim(twpe).xlim(xlim).A,1);
By = mean(pic.twpelim(twpe).xlim(xlim).By,1);
Bx = mean(pic.twpelim(twpe).xlim(xlim).Bx,1);



maxBy = max(By);

[I,J] = find(By>maxBy*0.5);


%plot(A,By,A(J),By(J),'.',[min(A) max(A)],0.5*maxBy*[1 1],'-')

plotyy(pic.zi,[Bx;By],pic.zi,A)
end


%% Check what attributes are missing for what time steps
% e.g.
missingAttr = no02m.get_missing_attributes('UK_ions');
% then use for it = missingAttr instead of it = 1:sim.length or it =
% 1:pic.nt
%% Energy partitioning, UB, UK, UT
pic = pic2;
times = pic.twci;
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB(it) = sum(Babs(:).^2)/2;  
  disp(sprintf('%g/%g',it,pic.nt))
end

%% Particle/plasma energy partitioning
sim = pic2;
tic;
clear UT UK
species_groups = {[1 3],[2 4]};
species_groups_str = {'ions','electrons'};
for it = 1:sim.length  
  sim_tmp = sim(it); 
  for iSpecies_group = 1:numel(species_groups)
    iSpecies = species_groups{iSpecies_group};
    species_str = species_groups_str{iSpecies_group};
    disp(sprintf('it = %g/%g, sp = %g ',it,sim.length,iSpecies))
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(iSpecies);
    %toc;
    pdyn = sim.mass(iSpecies(1))/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iSpecies_group) = 3/2*nansum(p(:));
    UK(it,iSpecies_group) = nansum(pdyn(:));
    %imagesc(sim.xi,sim.zi,squeeze(p)')
    %drawnow
    h5write_attr(sim_tmp,sim_tmp.twci,['UK_' species_str],UK(it,iSpecies_group))
    h5write_attr(sim_tmp,sim_tmp.twci,['UT_' species_str],UT(it,iSpecies_group))
  end
  toc
end