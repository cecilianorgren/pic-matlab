idist = 1;
timestep = 5000;
nss = 6;
for distnumber = 391:530%1:100;%300%101:200%1:100%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  distnumber
  read_sub_dir = '/1/';
  %txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
  txtfile = sprintf('/Users/cno062/tesla/cno062/df_cold_protons_n04/distributions/%05.0f/%.0f.dat',timestep,distnumber); % df04
  %txtfile = sprintf('/Users/cno062/tesla/cno062/df_cold_protons_n04/distributions/%.0f.dat',distnumber); % df04
  if not(exist(txtfile,'file'))
    error(sprintf('File not found: %s', txtfile))    
    %continue
  end  
    
  idist = idist + 1;
  %idist = distnumber;
  idist
  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile,nss);
    
    h = setup_subplots(2,nss);
    for isp = 1:nss
    hca = h(isp+nss);
    imagesc(hca,axes(:,isp),axes(:,isp),squeeze((fxy(:,:,isp)))');    
    hca = h(isp+0);
    imagesc(hca,axes(:,isp),axes(:,isp),squeeze((fxz(:,:,isp)))');
    end
    h(1).Title.String = sprintf('distnumber = %g',distnumber);
    h(2).Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlo,xhi,zlo,zhi);
    pause(0.1)
end
    