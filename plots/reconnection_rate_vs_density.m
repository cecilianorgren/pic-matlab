df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');

%% Calculate density at and above xline
df04n_Bxline_z1 = interp(df04n,df04n.x_xline,df04n.z_xline+1,df04n.twci,'Bx');
df04n_nxline_z1 = interp(df04n,df04n.x_xline,df04n.z_xline+1,df04n.twci,'ne');
df04n_nxline_z0 = interp(df04n,df04n.x_xline,df04n.z_xline+0,df04n.twci,'ne');

df04_Bxline_z1 = interp(df04,df04.x_xline,df04.z_xline+1,df04.twci,'Bx');
df04_nxline_z1 = interp(df04,df04.x_xline,df04.z_xline+1,df04.twci,'ne');
df04_nxline_z0 = interp(df04,df04.x_xline,df04.z_xline+0,df04.twci,'ne');

%% Figure
ylim_R = [0 0.15];
ytick_R = [0 0.05 0.1 0.15];

ylim_n = [0 1];

nrows = 2;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % ER, n
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,df04n.twci,([df04n_nxline_z0,df04n_nxline_z1])',df04n.twci,df04n.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end
if 1 % ER, B/sqrt(n)
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,df04n.twci,1./sqrt(df04n_nxline_z1),df04n.twci,df04n.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = [0 2];ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end

if 0 % Er, n
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,df04.twci,[df04_nxline_z0,df04_nxline_z1]',df04.twci,df04.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end

fig = gcf;
hlinks = linkprop(fig.Children,{'XLim'});
compact_panels(0.01)

for ip = 1:npanels
  h(ip).XGrid = 'on';
  %h(ip).YGrid = 'on';  
end