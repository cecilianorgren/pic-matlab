% Make example structure

% Radial coordinate system
theta = 1:2:360; % degrees
r = linspace(0,5,100);
[R,THETA] = meshgrid(r,theta);
nTheta = numel(theta);
nR = numel(r);

% Cartesian coordinate system
X = R.*cosd(THETA);
Y = R.*sind(THETA);

% Radial profile
fun_profile = @(r) r.*exp(-((r-2)/0.5).^2);
fun_profile_theta = @(theta) 1-0.2*cosd(theta);

radialprofile1D = fun_profile(r);
angularprofile = fun_profile_theta(theta);
profile_model = 1;
switch profile_model
  case 1 % simple smooth profile
    radialprofile2D = fun_profile(R);
  case 2 % simple radial profile with added noise
    radialprofile2D = fun_profile(R) + 0.5*rand(size(R)); % add some noise
  case 3 % add some theta dependence
    radialprofile2D = fun_profile(R).*fun_profile_theta(THETA); % add some noise    
  case 4 % add some theta dependence + noise
    radialprofile2D = fun_profile(R).*fun_profile_theta(THETA) + 0.5*rand(size(R));
  case 5 % add some theta dependence + noise
    radialprofile2D_1 = fun_profile(R).*fun_profile_theta(THETA) + 0.5*rand(size(R));  
    radialprofile2D_2 = fun_profile(R*1.5).*fun_profile_theta(THETA-180) + 0.5*rand(size(R));  
    radialprofile2D = radialprofile2D_1 + radialprofile2D_2;
end

nSavePeaks = 3;
colors = [     0    0.4470    0.7410;
          0.8500    0.3250    0.0980;
          0 0 0
          0.9290    0.6940    0.1250];
peak_location = nan(nTheta,nSavePeaks); % save 3 first peaks
% Find equation maximum value of the data
for iTheta = 1:nTheta % step through each angle theta
  radialprofile_tmp = radialprofile2D(iTheta,:);
  % Find peak using function findpeaks
  [PKS,LOCS] = findpeaks(radialprofile_tmp,r,'sort','descend');
  peak_location(iTheta,1:min(nSavePeaks,numel(LOCS))) = LOCS(1:min(nSavePeaks,numel(LOCS)));
end
% Cartesian coordinates
xpeak = peak_location.*cosd(theta');
ypeak = peak_location.*sind(theta');

% Plot
nrows = 3;
ncols = 1;
isub = 0;
for irows = 1:nrows
  for icols = 1:ncols
    isub = isub + 1;
    h(isub) = subplot(nrows,ncols,isub);
  end
end

% Plot
isub = 1;
if 0
  hca = h(isub); isub = isub + 1;  
  plot(hca,radius,radialprofile1D)
  hca.Title.String = 'Radial profile, without noise';
  hca.XLabel.String = 'r';
  %hca.XLabel.String = 'y';
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,X,Y,radialprofile2D)
  colormap('gray')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hca.Title.String = 'Data, three first peaks';
  hca.YLabel.String = 'x';
  hca.XLabel.String = 'y';
  % Add ocation of peaks
  hold(hca,'on')
  %plot(hca,xpeak,ypeak,'color',[0 0 0])
  %for iPeak = 1:nSavePeaks
  %  plot(hca,xpeak(:,iPeak),ypeak(:,iPeak),'linestyle','-','color',colors(iPeak,:))
  %end
  plot(hca,xpeak,ypeak,'linestyle','-','linewidth',1,'marker','.')
  hold(hca,'off')
end
if 1
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,X,Y,radialprofile2D)
  colormap('gray')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hca.Title.String = 'Data, main peaks';
  hca.YLabel.String = 'x';
  hca.XLabel.String = 'y';
  % Add ocation of peaks
  hold(hca,'on')
  %plot(hca,xpeak,ypeak,'color',[0 0 0])
  %for iPeak = 1:nSavePeaks
  %  plot(hca,xpeak(:,iPeak),ypeak(:,iPeak),'linestyle','-','color',colors(iPeak,:))
  %end
  plot(hca,xpeak(:,1),ypeak(:,1),'linestyle','-','linewidth',1,'marker','.')
  hold(hca,'off')
end
if 1 % peak location as a function of angle
  hca = h(isub); isub = isub + 1;  
  plot(hca,theta,peak_location)
  hca.Title.String = 'Location of peaks';
  hca.YLabel.String = 'r';
  hca.XLabel.String = 'theta';
end


%%
Bn_fun = @(lat,long) long;

long = linspace(0,360,100);
lat = linspace(0,90,50);
[LAT,LONG] = meshgrid(lat,long);
Bn = Bn_fun(LAT,LONG);

pcolor(lat,long,Bn)
shading flat

boundary_lat_fun = @(long) 70+10*cosd(long);

boundary_lat = boundary_lat_fun(long);
hold on;
plot(long,boundary_lat)

[BOUND,LONG] = meshgrid(boundary_lat,long);


