% Stream functions
c_eval('Se?.xz = vector_potential(x,z,ve?.x,ve?.z);',1:2) % stream function
c_eval('Si?.xz = vector_potential(x,z,vi?.x,vi?.z);',1:2) % stream function



nrows = 2;
ncols = 2;
npanels = nrows*ncols;
clear h;
isub = 1;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end

if 1 % vi1.x
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x,z,vi1.x')
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';  
end
if 1 % vi1.z
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x,z,vi1.z')
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';  
end
if 1 % ion stream function
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x,z,Si1.xz')
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';  
end