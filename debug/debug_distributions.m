% debug_distributions
%%
clear boxes
idist = 0;
for distnumber = 201%1:100%:10%:281%:20%40%10:281%5:280
  idist = idist + 1;
  %distnumber = 0;
  timestep = 8000;
  sub_dir = '/2/';


  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,sub_dir,distnumber); 
  %txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%.0f.dat',distnumber); % michael's perturbation  
  %txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,sub_dir,distnumber);
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%.0f.dat',distnumber); 
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = read_distributions(txtfile);
  vx = axes;
  vy = axes;
  vz = axes;
  imagesc(squeeze(vx(:,1)),squeeze(vz(:,1)),squeeze(fxz(:,:,1)))
  hca = gca;
  hca.Title.String = txtfile;
  hca.Title.Interpreter = 'none';
  colorbar
  pause(0.1)
  boxes(idist,:) = [xlo xhi zlo zhi];
end
box_size = (boxes(:,2)-boxes(:,1)).*(boxes(:,4)-boxes(:,3));

%%
nboxes = size(boxes,1);

hca = subplot(1,1,1);
imagesc(hca,x,z,ve1.x')
for ibox = 1:nboxes
  hold(hca,'on')
  hpatch = patch(hca,boxes(ibox,[1 1 2 2]),boxes(ibox,[3 4 4 3]),[0.8 0.8 0.8]);
  hpatch.FaceAlpha = 0.5;
  hold(hca,'off')
  pause(0.1)
end
