%% Load h5-objects
dist = PICDist('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/dists.h5');
pic = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');

%% Draw a map of distributions
iSpecies = [3 5];
iIter = 1;
zlim = [0 5];
xlim = [];

%% Draw distributions above 2D map of some parameter

twpe = 5000;
it = 1; % 5000
zlim = [0 5];
xlim = [170 200];
vlim = [-1.5 1.5];


twpe = 8000;
it = 2; % 8000
zlim = [0 5];
xlim = [150 200];
vlim = [-1.5 1.5];

pic = df04.twpelim(twpe).xlim([150 210]).zlim([-10 10]);
dst = dist.xlim(xlim).zlim(zlim);

iSpecies = [3 5];
iIter = 1;

%[n3,jx3,jy3,jz3,pxx3,pxy3,pxz3,pyy3,pyz3,pzz3] = pic.njp([3]);
%t3 = (pxx3+pyy3+pzz3)/3./n3; 
[n5,jx5,jy5,jz5,pxx5,pxy5,pxz5,pyy5,pyz5,pzz5] = pic.njp([5]);
t5 = (pxx5+pyy5+pzz5)/3./n5; 
varstr = 't5'; varlim = [0 0.15]; varcolormap = pic_colors('candy');
%varstr = 'pic.n([5])';
%varstr = 'pic.n([1 3 5])';
%varstr = 'pic.vx([5])'; varlim = [-1.2 1.2]; varcolormap = pic_colors('candy2');


nrows = 2;
ncols = 3;
h(1) = subplot(nrows,ncols,1:3);
for isub = 2:((nrows-1)*ncols+1) 
  h(isub) = subplot(nrows,ncols,(ncols-1)+isub);
  axis(h(isub),'square')
end

isub = 1;
if 1 % Colormap of some quantity
  hca = h(isub); isub = isub + 1;  
  var = eval(varstr);
  %imagesc(hca,pic.xi,pic.zi,pic.n([1 3 5])')
  imagesc(hca,pic.xi,pic.zi,var')  
  colormap(hca,varcolormap)
  hcb = colorbar('peer',hca);
  hca.YLabel.String = varstr;
  hold(hca,'on')
  iAx = 1:5:pic.nx;
  iAz = 1:5:pic.nz;
  Atmp = pic.A; 
  contour(hca,pic.xi(iAx),pic.zi(iAz),Atmp(iAx,iAz)',-25:1:0,'k')
  hold(hca,'off')
  hca.CLim = varlim;
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
end

for id = dst.indices{it}
  % Distributions
  isub = 2;
  if 1 % fxy     
    hca = h(isub); isub = isub + 1;
    f = dist.fxyz(it,id,iSpecies,3);
    imagesc(hca,f.v,f.v,f.f')
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hca.YDir = 'normal';   
    hca.XLabel.String = 'v_x (v_{A0})';
    hca.YLabel.String = 'v_y (v_{A0})';
  end
  if 1 % fxz
    hca = h(isub); isub = isub + 1;
    f = dist.fxyz(it,id,iSpecies,2);
    imagesc(hca,f.v,f.v,f.f')
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hca.YDir = 'normal';
    hca.XLabel.String = 'v_x (v_{A0})';
    hca.YLabel.String = 'v_z (v_{A0})';
  end
  if 1 % fyz
    hca = h(isub); isub = isub + 1;
    f = dist.fxyz(it,id,iSpecies,1);
    imagesc(hca,f.v,f.v,f.f')
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hca.YDir = 'normal'; 
    hca.XLabel.String = 'v_y (v_{A0})';
    hca.YLabel.String = 'v_z (v_{A0})';
  end
  hlinks = linkprop(h(2:4),'CLim','XLim','YLim');
  h(2).XLim = vlim;
  h(2).YLim = vlim;
  %h(2).CLim(1) = 0;
  %h(2).CLim(2) = 1e-2;  
  c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',2:4)
  
  % Plot location on map
  if exist('hbox','var')
    delete(hbox)
  end
  hold(h(1),'on')
  hbox = plot(h(1),f.x([1 2 2 1 1]),f.z([1 1 2 2 1]),'k');  
  
  hold(h(1),'off')
  pause(0.1)
end
