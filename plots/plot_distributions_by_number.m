for distnumber = 1:139%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  read_sub_dir = '/1/';
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
  if not(exist(txtfile,'file'))
    warning(sprintf('File not found: %s', txtfile))
    continue
  end  
    
  idist = idist + 1;
  %idist = distnumber;
  
  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile);
    
    h = setup_subplots(2,4);
    for isp = 1:4
    hca = h(isp+4);
    imagesc(hca,vx(:,isp),vy(:,isp),squeeze((fxy(:,:,isp)))');
    hca = h(isp+0);
    imagesc(hca,vx(:,isp),vz(:,isp),squeeze((fxz(:,:,isp)))');
    end
    pause(0.1)
end
    