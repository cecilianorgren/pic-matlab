J = 0.5;
n = 1;
vi = @(teti) (J/n)./(1+teti);
ve = @(teti) vi(teti).*(-teti);

tite = 1:0.5:10;
teti = sort(1./tite);

hca = subplot(1,1,1);
plot(hca,tite,[vi(1./tite);ve(1./tite)])
hca.YTick = -10:0.2:10;
grid on

