% use load_resaved_data.m
t = timesteps/50;
x = sim_info.x-mean(sim_info.x);
z = sim_info.z;
[Tx,X] = ndgrid(t,x);
[Tz,Z] = ndgrid(t,z);
% make_time_series

lA = -25:0;

%%
xpick = 0;
ix = find_closest_ind(x,xpick);
rx = 0;
ix = ix + [-rx:1:rx];
iz = 1:numel(z);

h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,t,z,squeeze(mean(A(:,ix,iz),2))');
shading(hca,'flat')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'z (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'A';

hold(hca,'on')
C = squeeze(mean(A(:,ix,iz),2));
contour(hca,Tz,Z,C,lA,'k')
hold(hca,'off')
hca.CLim = [-25 0];

colormap(hca,pic_colors('candy'))

hca = h(isub); isub = isub + 1;
icomp = 1;
imagesc(hca,t,z,squeeze(mean(B(:,ix,iz,icomp),2))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'z (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Bx';

hold(hca,'on')
C = squeeze(mean(A(:,ix,iz),2));
contour(hca,Tz,Z,C,lA,'k')
hold(hca,'off')

hca.CLim = 1*[-1 1];
colormap(hca,pic_colors('blue_red'))

hca = h(isub); isub = isub + 1;
icomp = 2;
imagesc(hca,t,z,squeeze(mean(E(:,ix,iz,icomp),2))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'z (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Ey';

hold(hca,'on')
C = squeeze(mean(A(:,ix,iz),2));
contour(hca,Tz,Z,C,lA,'k')
hold(hca,'off')

hca.CLim = 0.2*[-1 1];
colormap(hca,pic_colors('blue_red'))

%%
%zlim = [-10 10];
%ipx1 = find(x>xlim(1),1,'first');
%ipx2 = find(x<xlim(2),1,'last');
%ipz1 = find(z>zlim(1),1,'first');
%ipz2 = find(z<zlim(2),1,'last');

zpick = 0;
iz = find_closest_ind(z,zpick);
ix = 1:numel(x);

h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,t,x(ix),squeeze(A(:,ix,iz))')
shading(hca,'flat')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'x (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'A';

hold(hca,'on')
C = squeeze(A(:,ix,iz));
contour(hca,Tx,X,C,lA,'k')
hold(hca,'off')

colormap(hca,pic_colors('candy'))
hca.CLim = [-25 0];

hca = h(isub); isub = isub + 1;
icomp = 3;
imagesc(hca,t,x,squeeze(B(:,ix,iz,icomp))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'x (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Bz';

hold(hca,'on')
C = squeeze(A(:,ix,iz));
contour(hca,Tx,X,C,lA,'k')
hold(hca,'off')

hca.CLim = 0.5*[-1 1];
colormap(hca,pic_colors('blue_red'))

hca = h(isub); isub = isub + 1;
icomp = 2;
imagesc(hca,t,x,squeeze(E(:,ix,iz,icomp))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'x (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Ey';

hold(hca,'on')
C = squeeze(A(:,ix,iz));
contour(hca,Tx,X,C,lA,'k')
hold(hca,'off')

hca.CLim = 0.2*[-1 1];
colormap(hca,pic_colors('blue_red'))

hlink = linkprop(h,{'XLim','YLim'});
hlink.Targets(1).YLim = 60*[-1 1];
