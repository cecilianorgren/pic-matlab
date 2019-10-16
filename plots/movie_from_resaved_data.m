% first use load_resaved_data

nframes = size(E_ts,1);
cell_movie = cell(2,1);

x = sim_info.x-mean(sim_info.x);
z = sim_info.z;
xlim = [-30 30];
zlim = [-6 6]; 
ix1 = find(x>xlim(1),1,'first');
ix2 = find(x<xlim(2),1,'last');
iz1 = find(z>zlim(1),1,'first');
iz2 = find(z<zlim(2),1,'last');
ipx = ix1:1:ix2;
ipz = iz1:1:iz2;


doA = 1;
if doA    
  cA = 0*[0.8 0.8 0.8];
  cA = 0*[0.7 0.7 0.7];
  nA = 40;
  nA = [0:-0.5:min(A(:))];
  ipxA = ix1:5:ix2;
  ipzA = iz1:5:iz2;
end

h = setup_subplots(2,1);
clear hb
npanels = numel(h);
for iframe = 1:nframes
  % make figure
  isub = 1;
  if 1 % Ey
    hca = h(isub); isub = isub + 1; 
    variable = E_ts;
    %pcolor(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,3))'); 
    imagesc(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,2))'); 
    hcb = colorbar('peer',hca);
    hb(isub) = hcb; 
    hcb.YLabel.String = 'E.y';
    %shading(hca,'flat'); 
    colormap(hca,pic_colors('blue_red')); 
    hca.CLim = 0.5*[-1 1];
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    %hca.Title.String = sprintf('t\omega_{ce}=%0.5g',timesteps);
    hca.Title.String = sprintf('tw_{ce} = %05.0f',timesteps(iframe));
  end
  if 1 % vex
    hca = h(isub); isub = isub + 1; 
    variable = ve12_ts;
    %pcolor(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,3))'); 
    imagesc(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,1))'); 
    hcb = colorbar('peer',hca);
    hb(isub) = hcb; 
    hcb.YLabel.String = 've.x';
    %shading(hca,'flat'); 
    colormap(hca,pic_colors('blue_red')); 
    hca.CLim = 5*[-1 1];
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    %hca.Title.String = sprintf('t\omega_{ce}=%0.5g',timesteps);
    hca.Title.String = sprintf('tw_{ce} = %05.0f',timesteps(iframe));
  end
  if 0
    hca = h(isub); isub = isub + 1; 
    variable = nine;
    %pcolor(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz))'); 
    imagesc(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz))'); 
    %shading(hca,'flat'); 
    colormap(hca,pic_colors('blue_red')); 
    hca.CLim = 0.1*[-1 1];
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    hca.Title.String = 'ni-ne';
  end
  if 0
    hca = h(isub); isub = isub + 1;
    variable = je1+je2;
    pcolor(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,1))'); 
    shading(hca,'flat'); 
    colormap(hca,pic_colors('blue_red')); 
    hca.CLim = 0.5*[-1 1];
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    hca.Title.String = 'je.x';
  end
  if 0
    hca = h(isub); isub = isub + 1;
    variable = je1+je2;
    pcolor(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,3))'); 
    shading(hca,'flat'); 
    colormap(hca,pic_colors('blue_red')); 
    hca.CLim = 0.25*[-1 1];
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    hca.Title.String = 'je.z';
  end
  if 0
    hca = h(isub); isub = isub + 1;
    variable = ji1+ji2;
    pcolor(hca,x(ipx),z(ipz),squeeze(variable(iframe,ipx,ipz,3))'); 
    shading(hca,'flat'); 
    colormap(hca,pic_colors('blue_red')); 
    hca.CLim = 0.25*[-1 1];
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    hca.Title.String = 'ji.z';
  end
  
  for ipanel = 1:npanels
    h(ipanel).Box = 'on';
    %colorbar('peer',h(ipanel));
    if doA
      hold(h(ipanel),'on')
      hcont = contour(h(ipanel),x(ipxA),z(ipzA),squeeze(A_ts(iframe,ipxA,ipzA))',nA,'color',cA,'linewidth',0.7);  
      hold(h(ipanel),'off')  
    end
  end
  
  % collect frames
  currentBackgroundColor = get(gcf,'color');
  set(gcf,'color',[1 1 1]);
  drawnow      
  tmp_frame = getframe(gcf);
  %cell_movies{imovie}(itime) = tmp_frame;
  if iframe == 1 % initialize animated gif matrix
    [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
    %map(end+1,:) = get(gcf,'color');
    im_tmp(1,1,1,nframes) = 0;
    cell_movie{1} = map;
    cell_movie{2} = im_tmp;
  else
    cell_movie{2}(:,:,1,iframe) = rgb2ind(tmp_frame.cdata,cell_movie{1},'nodither');
  end       
      pause(0.1)
end

% make gif
% imwrite(cell_movie{2},cell_movie{1},[savedir_root 'Ey_vex.gif'],'DelayTime',0.0,'LoopCount',0)