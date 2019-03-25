% Physical constants
%units = irf_units;
me = 9.1094e-31;
mp = 1.6726e-27;
mass = me;
e = 1.6022e-19;

% distribution integrated over two (perpendicular) directions already
% Tpar = Tperp1 = Tperp2 assumed
f_red = @(v,n,T,vd,m) n/((pi)^(1/2)*sqrt(2*e*T/m))*exp(-(v-vd).^2./(2*e*T/m));
f_T = @(v,m) m.*v.^2/(2*e);
f_vt = @(T,m) sqrt(2*e*T./m);

% Distribution properties
iexample = 1;
switch iexample
  case 1 % two distributions, same temperature, one still one drifting  
    m = [1 1]*mass; % kg
    n = [1 1]*1e-6; n_tot = sum(n); % m^-3
    % vt from T
    T = [10 10]; % eV    
    vt = f_vt(T,m); % m/s
    % T from vt
    vt = 0.5*[1000 1000]*1e3; % m/s
    T = f_T(vt,m); % eV
    vd = [0 10000]*1e3; % m/s
    ns = numel(n); % number of populations
  case 2 % two distributions, same temperature, opposite drift speed
    m = [1 1]*mass; % kg
    n = [1 1]*1e-6; n_tot = sum(n); % m^-3
    % vt from T
    T = [10 10]; % eV    
    vt = f_vt(T,m); % m/s
    % T from vt
    vt = 0.5*[1000 1000]*1e3; % m/s
    T = f_T(vt,m); % eV
    vd = [-5000 5000]*1e3; % m/s
    ns = numel(n); % number of populations
  case 3 %
    m = [1 1 1]*mass; % kg
    n = [1 1 0.2]*1e-6; n_tot = sum(n); % m^-3
    % vt from T
    T = [10 10 10]; % eV    
    vt = f_vt(T,m); % m/s
    % T from vt
    vt = [10000 1000 500]*1e3; % m/s
    T = f_T(vt,m); % eV
    vd = [-5000 5000 8000]*1e3; % m/s
    ns = numel(n); % number of populations
  case 4 % bug check
    m = [1 1]*mass; % kg
    n = [2 1]*1e-6; n_tot = sum(n); % m^-3
    % vt from T
    T = [10 10]; % eV    
    vt = f_vt(T,m); % m/s
    % T from vt
    vt = [500 500]*1e3; % m/s
    T = f_T(vt,m); % eV
    vd = [0 000]*1e3; % m/s
    ns = numel(n); % number of populations
end
% Constructing f
nv = 2000;
v = linspace(-30000,30000,nv)*1e3; % if the distributions go outside the 
                                   % defined v-interval, the recalculated 
                                   % moments will be less accurate
dv = v(2)-v(1);

legends_f = {};
f_sum_str = 'fsum =';
f_sep = zeros(ns,nv); % different populations separate
f_sep_vph = zeros(ns,1); 

for is = 1:ns
  f_sep(is,:) = f_red(v,n(is),T(is),vd(is),m(is));
  legends_f{is} = sprintf('f%g: n = %g cm-3, T = %g eV, vt = %g km/s, vd = %g km/s',is,n(is)*1e6,T(is),vt(is)*1e-3,vd(is)*1e-3);
  if is == 1, f_sum_str = sprintf('%s f%g',f_sum_str,is);
  else, f_sum_str = sprintf('%s + f%g',f_sum_str,is); end
end
f_sum = sum(f_sep,1); % sum of populations

% integrate f to get common moments
n_tot = sum(f_sum)*dv; % m^-3
vd_tot = sum(f_sum.*v)*dv./n_tot; % m/s
if 0
  vt_tot = sqrt(sum(abs(v-vd_tot).^2.*f_sum)*dv*2/n_tot); % m/s
else
  % n^-1*int((v-vd)^2*f))=vt^2/2
  % vt^2 = sqrt(2/m)*int((v-vd)^2*f))
  vt2_tot = (2/n_tot)*sum(abs(v-vd_tot).^2.*f_sum)*dv; % (m/s)^2
  vt_tot = sqrt(vt2_tot); % m/s
end
T_tot = mass*vt_tot^2/(2*e); % eV
f_tot = f_red(v,n_tot,T_tot,vd_tot,mass);

% Plotting
linewidth = 1;
%color_tr = [1 0.9 0.4];
color_sum = [0 0 0];
color_tot = [0.7 0.7 0.7];
color_sep = [0.0000    0.4470    0.7410;...
             0.8500    0.3250    0.0980;...
             0.9290    0.6940    0.1250;...
             0.4940    0.1840    0.5560;...
             0.4660    0.6740    0.1880;...
             0.3010    0.7450    0.9330;...
             0.6350    0.0780    0.1840];
linestyle_sum = '--';
linestyle_tot = '-';
linestyle_sep = '-';

hca = subplot(1,1,1);
% plot f_sep
set(hca,'ColorOrder',color_sep,'ColorOrderIndex',1)
h_pl_sep = plot(hca,v*1e-3,f_sep,'linewidth',linewidth,'linestyle',linestyle_sep);

hold(hca,'on')
% plot f_sum
h_pl_tot = plot(hca,v*1e-3,f_sum,'linewidth',linewidth,'linestyle',linestyle_sum ,'color',color_sum);
legends_f{end+1} = f_sum_str;

% plot f_tot, recalculated maxwellian from f_sum
h_pl_tot = plot(hca,v*1e-3,f_tot,'linewidth',linewidth,'linestyle',linestyle_tot,'color',color_tot);
legends_f{end+1} = sprintf('ftot: n = %g cm-3, T = %g eV, vt = %g km/s, vd = %g km/s',n_tot*1e6,T_tot,vt_tot*1e-3,vd_tot*1e-3);

hold(hca,'off')

legend(hca,legends_f{1:(ns+2)},'location','northoutside')
hca.XLabel.String = 'v (km/s)';
hca.YLabel.String = {'f (m-4s-1)','(reduced to 1D)'};
