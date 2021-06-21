%% Load data
twpe = 24000; % timestep in units of inverse electron plasma frequencies, 
              % twci = twpe/(wpewce*mime), wpewce = 2, mime = 100
              % twpe = 18000, twci = 090
              % twpe = 20000, twci = 100
              % twpe = 22000, twci = 110
              % twpe = 24000, twci = 120
              % twpe = 25000, twci = 115
filepath = sprintf('/Path/to/your/file/fields-%5.0f.dat',twpe);
data = function_load_data(filepath); % load data into struct (for help type
                                     % 'help struct' in terminal
%% Access data, examples
n_hot = data.hot_ion.n;            % density for hot Harris sheet population
n_cold = data.cold_ion.n;          % density for cold inflow populations 
n_cold_top = data.cold_ion_top.n;  % density for cold inflow population, originating from the north
n_cold_bot = data.cold_ion_bot.n;  % density for cold inflow population, originating from the south

%% Plot, example 
nrows = 7;
ncols = 1;
ipanel = 1;

% Magnetic field component Bz
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
imagesc(hca,data.x_di,data.z_di,data.Bz');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'B_z';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';

% Electric field component Ey
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
imagesc(hca,data.x_di,data.z_di,data.Ey');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'E_y';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';

% Current density Jy
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
Jy = data.hot_ion.jy - data.hot_ele.jy + data.cold_ion.jy - data.cold_ele.jy;
imagesc(hca,data.x_di,data.z_di,Jy');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'J_y';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';
hca.CLim = [-1 1];

% Density of cold ions originating from the north
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
imagesc(hca,data.x_di,data.z_di,data.cold_ion_top.n');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'n_{i,cold,top}';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';

% Pressure of cold ions originating from the north
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
imagesc(hca,data.x_di,data.z_di,data.cold_ion_top.p');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'p_{i,cold,top}';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';

% Pressure of cold ions originating from north and south
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
imagesc(hca,data.x_di,data.z_di,data.cold_ion.p');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'p_{i,cold}';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';

% Velocity component vx of cold ions originating from north and south
hca = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
imagesc(hca,data.x_di,data.z_di,data.cold_ion.vx');
hca.YDir = 'normal';
hcb = colorbar('peer',hca); hcb.YLabel.String = 'v_{x,i,cold}';
hca.XLabel.String = 'x/d_i';
hca.YLabel.String = 'z/d_i';
hca.CLim = [-1.5 1.5];