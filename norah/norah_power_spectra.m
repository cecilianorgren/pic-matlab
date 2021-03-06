%% Load data
filename = '/Users/cecilia/MATLAB/pic-matlab/norah/Haloween2003.mat';
load(filename)

% See content
s = whos('-file',filename);
varnames = {s.name};

% Make time into EpochTT object (that's what irfu-matlab uses)
% ISO format, 2018-07-05T20:19:00.00Z
nt = numel(Year);
[Year' Month' Day' Hour' Min' Sec'];
time = [];
for it = 1:nt
  iso_str_time = sprintf('%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0fZ',Year(it),Month(it),Day(it),Hour(it),Min(it),Sec(it));  
  time = [time EpochTT(iso_str_time)];
end

% Make TSeries objects (that's what irfu-matlab uses)
tBd_obs = irf.ts_scalar(time,Obs_dbd);
tBe_obs = irf.ts_scalar(time,Obs_dbe);
tBn_obs = irf.ts_scalar(time,Obs_dbn);

tBd_rusbx = irf.ts_scalar(time,Model_rusbx_dbd);
tBe_rusbx = irf.ts_scalar(time,Model_rusbx_dbe);
tBn_rusbx = irf.ts_scalar(time,Model_rusbx_dbn);

tBd_sokbx = irf.ts_scalar(time,Model_sokbx_dbd);
tBe_sokbx = irf.ts_scalar(time,Model_sokbx_dbe);
tBn_sokbx = irf.ts_scalar(time,Model_sokbx_dbn);

tBd_sokbx0 = irf.ts_scalar(time,Model_sokbx0_dbd);
tBe_sokbx0 = irf.ts_scalar(time,Model_sokbx0_dbe);
tBn_sokbx0 = irf.ts_scalar(time,Model_sokbx0_dbn);

if 1 % Test plot
  h = irf_plot(1);
  if 1
    hca = irf_panel('dBd');
    irf_plot(hca,{tBd_obs,tBd_rusbx,tBd_sokbx,tBd_sokbx0},'comp')
    hca.YLabel.String = 'dBd';
    irf_legend(hca,{'Obs','Rus','Sok','Sok Bx=0'},[0.98 0.98])
  end
end

%% Make spectrogram
dt = time(2)-time(1);
samplFreqHz = 1/dt;
nFft = tBd_obs.length; % 1 min?
sBd_obs = irf_powerfft(tBd_obs,nFft,samplFreqHz);
sBd_rusbx = irf_powerfft(tBd_rusbx,nFft,samplFreqHz);
sBd_sokbx = irf_powerfft(tBd_sokbx,nFft,samplFreqHz);
sBd_sokbx0 = irf_powerfft(tBd_sokbx0,nFft,samplFreqHz);
varstrs = {'sBd_obs','sBd_rusbx','sBd_sokbx','sBd_sokbx0'};

% Smooth the data in logarithimcally spaced bins so that the log spectra is
% easer to read.
% units: B^2/Hz, to get the total power in a bin, need to multiply with bin
% width, and after collecting bins, need to again divide by the new bin
% width. HAD NOT IMPLEMENTED THIS BEFORE.

% original frequency bins
f_orig = sBd_sokbx0.f;
df_orig = diff(f_orig(1:2)); % uniform
% f_orig_edges = [f_orig(1)-0.5*df_orig(1); f_orig(1:end-1)+0.5*df_orig; f_orig(end)+0.5*df_orig(end)];
% if f_orig_edges(1)<0, f_orig_edges(1) = 0; end
% new frequency bins
nf = 100;
fsmooth_edges = logspace(-5,-2,nf+1);
fsmooth_edges = linspace(1e-5,1e-2,nf+1);
df = diff(fsmooth_edges);
fsmooth = fsmooth_edges(1:end-1)+0.5*df;
psmooth = zeros(1,nf);
for ivar = 1:numel(varstrs)
  var = eval(varstrs{ivar});
  bins = discretize(var.f,fsmooth_edges);
  for ibin = 1:nf
    ptmp = sum(var.p{1}(find(bins==ibin)))*df_orig;
    psmooth(ibin) = ptmp/df(ibin);
  end
  var.fsmooth = fsmooth;
  var.psmooth{1} = psmooth;
  % make step-like resampled spectra, to show bin edges
  var.fstep = reshape([fsmooth_edges(1:end-1);fsmooth_edges(2:end)],(numel(fsmooth_edges)-1)*2,1)';
  var.pstep{1} = reshape([psmooth;psmooth],numel(psmooth)*2,1)';
  
  eval([varstrs{ivar} '= var;'])
end

if 1 % Test plot
  nrows = 3;
  ncols = 1;
  npanels = nrows*ncols;
  h = setup_subplots(nrows,ncols);
  isub = 1;
  if 1
    hca = h(isub); isub = isub + 1;
    loglog(hca,sBd_obs.f,sBd_obs.p{1},...
               sBd_rusbx.f,sBd_rusbx.p{1},...
               sBd_sokbx.f,sBd_sokbx.p{1},...
               sBd_sokbx0.f,sBd_sokbx0.p{1}...
      )
    hca.XLabel.String = 'Frequency';
    hca.YLabel.String = 'Power';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    irf_legend(hca,{'Obs','Rus','Sok','Sok Bx=0'}',[0.02 0.02])
  end
  if 1 % smoothed data, compare single one
    hca = h(isub); isub = isub + 1;
    loglog(hca,sBd_obs.f,sBd_obs.p{1},...
               sBd_obs.fsmooth,sBd_obs.psmooth{1},'.-',...
               sBd_obs.fstep,sBd_obs.pstep{1},'-'...
      )
    hca.XLabel.String = 'Frequency';
    hca.YLabel.String = 'Power Obs';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    irf_legend(hca,{'original','smoothed (bin center)','smoothed (bin width)'}',[0.98 0.98])
  end
  if 1 % smoothed data, all
    hca = h(isub); isub = isub + 1;
    loglog(hca,sBd_obs.fsmooth,sBd_obs.psmooth{1},'.-',...
               sBd_rusbx.fsmooth,sBd_rusbx.psmooth{1},'.-',...
               sBd_sokbx.fsmooth,sBd_sokbx.psmooth{1},'.-',...
               sBd_sokbx0.fsmooth,sBd_sokbx0.psmooth{1},'.-'...
      )
    hca.XLabel.String = 'Frequency';
    hca.YLabel.String = 'Power';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    irf_legend(hca,{'Obs','Rus','Sok','Sok Bx=0'}',[0.02 0.02])
  end
end

%% Filter data and compare
fs = samplFreqHz;
flow = 0;
fhigh = 5e-4;
forder = 5;
if 1 % Test plot
  h = irf_plot(2);
  if 1
    hca = irf_panel('dBd');
    irf_plot(hca,{tBd_obs,...
                  tBd_rusbx,...
                  tBd_sokbx,...
                  tBd_sokbx0},'comp')
    hca.YLabel.String = 'dBd';
    irf_legend(hca,{'Obs','Rus','Sok','Sok Bx=0'},[0.98 0.98])
  end
  if 1
    hca = irf_panel('dBd filt');    
    irf_plot(hca,{tBd_obs.filt(flow,fhigh,fs,forder),...
                  tBd_rusbx.filt(flow,fhigh,fs,forder),...
                  tBd_sokbx.filt(flow,fhigh,fs,forder),...
                  tBd_sokbx0.filt(flow,fhigh,fs,forder)},'comp')
    hca.YLabel.String = {'dBd filt',sprintf('f = [%g,%g]',flow,fhigh)};
    irf_legend(hca,{'Obs','Rus','Sok','Sok Bx=0'},[0.98 0.98])
  end
end






