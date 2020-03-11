function [v,f] = construct_maxwellian_from_f(f)

doPlot = 1;

[fn,fv,fT] = calculate_moments(f); % Moments that go into a maxwellian distribution
% obs this temperature above is not scaled to mass, so it's bascially T/m
vt = 2*sqrt(fT);

fun_max = @(v,n,vd,T) n*(pi*vt.^3).^(-0.5).*exp(-0.5*(v-vd).^2./(vt.^2));

[VX,VY,VZ] = ndgrid(f.v,f.v,f.v);


if doPlot
  figure(77)  
  h = setup_subplots(1,1);
  isub = 1;
  
  if 1 % 
    hca = h(isub); isub = isub + 1;
    %%
    fplot = f.f;
    fplot(fplot==0) = NaN;
    h_scat = scatter3(hca,VX(:),VY(:),VZ(:),fplot(:)*0+5,fplot(:));
    h_scat.MarkerFaceColor = h_scat.MarkerEdgeColor;
    h_scat.MarkerEdgeColor = 'none';
    h_scat.MarkerFaceAlpha = 0.5;
    colormap(hca,pic_colors('candy'))
  end
  
  hca = h(isub); isub = isub + 1;
  
end