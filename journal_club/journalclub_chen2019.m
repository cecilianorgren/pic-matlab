L1 = [0.9311,-0.2407,0.2741];
M1 = [0.0806,0.8687,0.4889];
N1 = [-0.3558,-0.4330,0.8282];
L2 = [0.9910,-0.0513,-0.1240];
M2 = [0.0331,0.9889,-0.1445];
N2 = [0.1300,-0.1391,0.9817];
L3 = [0.9840,-0.1685,-0.0578];
M3 = [-0.0030,0.3090,-0.9511];
N3 = [0.1781,0.9360,0.3036];

colors = pic_colors('matlab');

hca = subplot(1,1,1);

quiver3(hca,0,0,0,L1(1),L1(2),L1(3),'color',colors(1,:),'linewidth',1,'displayname','L_1'); hold(hca,'on')
quiver3(hca,0,0,0,M1(1),M1(2),M1(3),'color',colors(1,:),'linestyle','--','displayname','M_1');
quiver3(hca,0,0,0,N1(1),N1(2),N1(3),'color',colors(1,:),'displayname','N_1');

quiver3(hca,0,0,0,L2(1),L2(2),L2(3),'color',colors(2,:),'linewidth',1,'displayname','L_2');
quiver3(hca,0,0,0,M2(1),M2(2),M2(3),'color',colors(2,:),'linestyle','--','displayname','M_2');
quiver3(hca,0,0,0,N2(1),N2(2),N2(3),'color',colors(2,:),'displayname','N_2');

quiver3(hca,0,0,0,L3(1),L3(2),L3(3),'color',colors(3,:),'linewidth',1,'displayname','L_3');
quiver3(hca,0,0,0,M3(1),M3(2),M3(3),'color',colors(3,:),'linestyle','--','displayname','M_3');
quiver3(hca,0,0,0,N3(1),N3(2),N3(3),'color',colors(3,:),'displayname','N_3'); hold(hca,'off')

axis(hca,'equal')
axis(hca,'square')
legend(hca)
hca.FontSize = 14;

hca.XLabel.String = 'x (GSM)';
hca.YLabel.String = 'y (GSM)';
hca.ZLabel.String = 'z (GSM)';

