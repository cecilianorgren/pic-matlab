J = 0.5;
n = 1;
vi = @(teti) (J/n)./(1+teti);
ve = @(teti) vi(teti).*(-teti);

tite = 0.001:0.5:10;
tite = logspace(-3,log10(10),100);
teti = sort(1./tite);

hca = subplot(1,1,1);
plot(hca,tite,[vi(1./tite);ve(1./tite)],'linewidth',1)
hca.YTick = -10:0.2:10;
grid on
hca.XLabel.String = 'T_i/T_e';
hca.YLabel.String = 'v_y';

hold(hca,'on')
%plot(hca,5,[vi(1./5);ve(1./5)],'*','linewidth',1)
plot(hca,5,[vi(1./5);ve(1./5)],'*k','linewidth',1)
hold(hca,'off')

legend(hca,{'v_i','v_e'},'box','off')


hca.LineWidth = 1;
hca.FontSize = 18;


