

h5filepath = '/Users/cecilia/Data/PIC/rec_onset_4/data_h5/fields.h5';
datapath = '/Users/cecilia/Data/PIC/rec_onset_4/data/';
%h5filepath = '/Users/cecilia/Data/PIC/rec_onset_asym/data_h5/fields.h5';
%datapath = '/Users/cecilia/Data/PIC/rec_onset_asym/data/';

h5filepath = '/Users/cecilia/Data/PIC/varying_inflow_density_steps_1/data_h5/fields.h5';
datapath = '/Users/cecilia/Data/PIC/varying_inflow_density_steps_1/data/';



h5filepath = '/Users/cno062/Data/PIC/rec_onset_4/data_F025_E005_TITE05/fields.h5';
datapath = '/Users/cno062/Data/PIC/rec_onset_4/data_F025_E005_TITE05/';

h5filepath = '/Users/cno062/Data/PIC/rec_onset_4/data_F025_E005_TITE10/fields.h5';
datapath = '/Users/cno062/Data/PIC/rec_onset_4/data_F025_E005_TITE10/';

h5filepath = '/Users/cno062/Data/PIC/varying_tite/tite_05/fields.h5';
datapath = '/Users/cno062/Data/PIC/varying_tite/tite_05/data/';

%h5filepath = '/Users/cno062/Data/PIC/varying_tite/tite_10/fields.h5';
%datapath = '/Users/cno062/Data/PIC/varying_tite/tite_10/data/';

nSpecies = 4;
h5write_fields(datapath,h5filepath,[0000:1000:20000],nSpecies)


%h5filepath = '/Users/cecilia/Data/PIC/rec_onset_4/data_F0_50_E005_TITE10/fields.h5';
%datapath = '/Users/cecilia/Data/PIC/rec_onset_4/data_F0_50_E005_TITE10/';

%nSpecies = 4;
%h5write_fields(datapath,h5filepath,[0000:1000:7000],nSpecies)



pic = PIC(h5filepath); % If you have many times saved, this can take up to a minute
h5write_fields_complement(pic)
  
pic = PIC(h5filepath);


%%
h5filepath = '/Users/Cecilia/Data/PIC/data/data_h5/fields.h5';
datapath = '/Users/Cecilia/Data/PIC/data/';
nSpecies = 4;
h5write_fields(datapath,h5filepath,600,nSpecies)
pic = PIC(h5filepath);

%% debug
txtfile1 = '/Users/cecilia/Data/PIC/rec_onset_3/data/fields-00001.dat';
%txtfile2 = '/Users/cecilia/Data/PIC/rec_onset_3/data/fields-00040.dat';
tic; [varstrs1,vars1] = read_data_no_normalization(txtfile1,nSpecies); toc  
%tic; [varstrs2,vars2] = read_data_no_normalization(txtfile2,nSpecies); toc

%%
nSpecies = 4;
for ifile = 1:40
  txtfile = sprintf('/Users/cecilia/Data/PIC/rec_onset_3/data/fields-000%02.0f.dat',ifile);
  [varstrs,vars] = read_data_no_normalization(txtfile,nSpecies);
  Ey = vars{find(contains(varstrs,'bz'))};
  imagesc(Ey')
  set(gca,'clim',1*[-1 1])
  title(gca,sprintf('ifile = %g',ifile))
  pause(0.1)
  
end

%% Check total driving
fontsize = 14;

pic_top_center = pic.xlim(mean(pic.xi)+0.1*[-1 1]).zlim(pic.zi(end)+[-0.01 0]);
Ey_at_boundary = squeeze(mean(mean(pic_top_center.Ey,1),2));
Flux_added = cumtrapz(pic.twci,Ey_at_boundary);
pic_inflow_center_t0 = pic(1).xlim(mean(pic.xi)+0.1*[-1 1]).zlim([0 pic.zi(end)]);
Bx_preexisting = squeeze(mean(pic_inflow_center_t0.Bx,1));
Flux_preexisting = trapz(pic_inflow_center_t0.zi,Bx_preexisting);

h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,pic.twci,Ey_at_boundary)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = sprintf('E_y^{drive} (B_0v_A)');
irf_legend(hca,sprintf('x = [%.2f,%.2f], z = [%.2f,%.2f]',pic_top_center.xi(1),pic_top_center.xi(end),pic_top_center.zi(1),pic_top_center.zi(end)),[0.98 0.98],'color','k')

if 0
plot(hca,pic.twci,Bx_preexisting)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'B_x';
irf_legend(hca,sprintf('x = [%.2f,%.2f], z = [%.2f,%.2f]',pic_inflow_center_t0.xi(1),pic_inflow_center_t0.xi(end),pic_inflow_center_t0.zi(1),pic_inflow_center_t0.zi(end)),[0.98 0.98],'color','k')
end  

hca = h(isub); isub = isub + 1;
plot(hca,pic.twci,Flux_added,pic.twci,pic.twci*0+Flux_preexisting,'--')
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = '\Phi (B_0d_i)';
hca.YLim(2) = ceil(Flux_preexisting);
irf_legend(hca,{'\Phi^{drive} = \int E_y^{drive}dt','\Phi^{pre-existing} = \int B_x(t=0) dz'}',[0.02 0.7])

hca = h(isub); isub = isub + 1;
plot(hca,pic.twci,Flux_added/Flux_preexisting)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = '\Phi^{drive}/\Phi^{pre-existing}';
%irf_legend(hca,{'\int E_y^{drive}dt','\int B_x(t=0) dz'}',[0.02 0.7])

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

compact_panels(0.02,0.00)

%% Estimating ts from E^peak
mime = 100;
wpewce = 2;
rel_added_flux = [0.25:0.25:1];
rel_added_flux = [0.1:0.1:5];
rel_added_flux = logspace(-2,2,10);
leg_str = arrayfun(@(x) ['\Phi_d/\Phi_0 = ' sprintf('%.2f',x)],rel_added_flux','UniformOutput',false);
B0 = 1;
L = 2;
zmax = 12.8;
%Ey_peak = 0.2;

ts = @(Ey_peak,rel_added_flux) rel_added_flux*B0*L*log(cosh(zmax/L))./(Ey_peak*exp(1));
%ts = @(Ey_peak,rel_added_flux) rel_added_flux*B0*zmax./(Ey_peak*exp(-1));
Ey_peak = linspace(0.01,10,200);

h = setup_subplots(1,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,Ey_peak,ts(Ey_peak,rel_added_flux(1)))
hold(hca,'on')
for ii = 2:numel(rel_added_flux)
  plot(hca,Ey_peak,ts(Ey_peak,rel_added_flux(ii)))
end
hold(hca,'off')
hca.XLabel.String = 'E_y^{peak} (B_0v_A)';
hca.YLabel.String = 't_s\omega_{ci}';
if 0
hca.YTick = 0:20:200;
hca.XTick = 0:0.1:10;
hca.YLim(2) = 100;
end
hleg = legend(hca,leg_str,'location','eastoutside');
hleg.Title.String = 'Added flux';


drawnow
c_eval('h(?).Position(2) = h(?).Position(2) + 0.08;',1:numel(h))
c_eval('h(?).Position(4) = h(?).Position(4) - 0.08;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))

c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

hca.YScale = 'log';
hca.XScale = 'log';

if 0
ax1 = hca;
ylim1 = ax1.YLim;
ytick1 = ax1.YTick;

yyaxis right
ax2 = gca;
ax2.YLim = ylim1*mime*wpewce;
ax2.YTick = ytick1*mime*wpewce;
ax2.YLabel.String = 't_s\omega_{pe} = (m_i/m_e)(\omega_{pe}/\omega_{ce})t_s\omega_{ci}';
drawnow
end