[saddle_locations,saddle_values] = saddle(A,'sort');
Ax = saddle_values(1);

if not(isfield(B,'abs')), B.abs = sqrt(B.x.^2 + B.y.^2 + B.z.^2); end
pe12.scalar = (pe12.xx+pe12.yy+pe12.zz)/3;
pi12.scalar = (pi12.xx+pi12.yy+pi12.zz)/3;


% Energy densities
UB = 0.5*B.abs.^2;
Ute = 3/2*pe12.scalar;
Uti = 3/2*pi12.scalar;
Uke = mass(2)/mass(1)*0.5*ne12.*(ve12.x.^2 + ve12.y.^2 + ve12.z.^2);
Uki = mass(1)/mass(1)*0.5*ni12.*(vi12.x.^2 + vi12.y.^2 + vi12.z.^2);


%%

vars = {UB,Uti,Ute,Uki,Uke};
nVars = numel(vars);
minA = -27; maxA = -15; nA = 100;
edgesA = linspace(minA,maxA,nA);
dA = edgesA(2) - edgesA(1);
levelsA = edgesA(2:end)+0.5*dA;
tic;
[volA,varsA] = calc_fluxtube_content(x,z,A,edgesA,vars);
toc;
allVars = nan(nA-1,nVars);
for iVar = 1:nVars
  allVars(:,iVar) = varsA{iVar};
end

%%
h = setup_subplots(3,2);
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  ipx = 1:5:numel(x);
  ipz = 1:5:numel(z);
  imagesc(hca,x(ipx),z(ipz),UB(ipx,ipz)')
  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U_B'; 
  hold(hca,'on')
  [hh,cc] = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',floor(edgesA(1)):1:ceil(edgesA(end)),'k');
  %clabel(hh,cc);
  hold(hca,'off')
  hca.CLim = [min(UB(:)),max(UB(:))];
end
if 1
  hca = h(isub); isub = isub + 1;  
  plot(hca,levelsA,volA,'linewidth',1.5)
  
  hca.XLabel.String = 'A';
  hca.YLabel.String = 'Fluxtube volume';
  hold(hca,'on')
  plot(hca,Ax*[1 1],hca.YLim,'k--')
  hold(hca,'off')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;
  %plot(hca,levelsA,varsA{1},levelsA,varsA{2},levelsA,varsA{3},levelsA,varsA{4},levelsA,varsA{5})
  plot(hca,levelsA,allVars,'linewidth',1.5)
  %plot(hca,levelsA,cumsum(allVars,2))
  hca.XLabel.String = 'A';
  hca.YLabel.String = 'Energy density';
  hold(hca,'on')
  plot(hca,Ax*[1 1],hca.YLim,'k--')
  hold(hca,'off')
  legend(hca,{'UB','Uti','Ute','Uki','Uke','A_{X line}'})
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;
  %plot(hca,levelsA,varsA{1},levelsA,varsA{2},levelsA,varsA{3},levelsA,varsA{4},levelsA,varsA{5})
  plot(hca,levelsA,allVars,'linewidth',1.5)
  %plot(hca,levelsA,cumsum(allVars,2))
  hca.XLabel.String = 'A';
  hca.YLabel.String = 'Energy density';
  hold(hca,'on')
  plot(hca,Ax*[1 1],hca.YLim,'k--')
  hold(hca,'off')
  legend(hca,{'UB','Uti','Ute','Uki','Uke','A_{X line}'})
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;
  %plot(hca,levelsA,varsA{1},levelsA,varsA{2},levelsA,varsA{3},levelsA,varsA{4},levelsA,varsA{5})
  %plot(hca,levelsA,allVars,'linewidth',1.5)
  plot(hca,levelsA,cumsum(allVars,2),'linewidth',1.5)
  hca.XLabel.String = 'A';
  hca.YLabel.String = 'Energy density';
  hold(hca,'on')
  plot(hca,Ax*[1 1],hca.YLim,'k--')
  hold(hca,'off')
  legend(hca,{'UB','Uti','Ute','Uki','Uke','A_{X line}'})
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;  
  plot(hca,levelsA,varsA{2}*teti,levelsA,varsA{3},'linewidth',1.5)
  %plot(hca,levelsA,cumsum(allVars,2))
  hca.XLabel.String = 'A';
  hca.YLabel.String = 'Energy density';
  hold(hca,'on')
  plot(hca,Ax*[1 1],hca.YLim,'k--')
  hold(hca,'off')
  legend(hca,{sprintf('Uti x %g',teti),'Ute','A_{X line}'})
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end