function varargout = read_data_no_normalization(txtfile,nss)  
% Reads data an returns two cell arrays, one with char variable
% name and one with data
% [varstrs,vars] = read_data_no_normalization(txtfile,nss);
%   nss - number of species
% Examples:
%   txtfile = '/Users/cecilia/tesla/cno062/guide_field_05/data/fields-00030.dat';
%   nss = 6;
%   [varstrs,vars] = read_data_no_normalization(txtfile,nss);


%% Load data
[fid, message] = fopen(txtfile,'r','ieee-le');
if fid < 0
  error('Failed to open file "%s" because "%s"', txtfile, message);
end
header = fread(fid,1,'integer*8');

it = fread(fid,1,'integer*4');                                % it
dt = fread(fid,1,'real*4');                                   % dt
teti = fread(fid,1,'real*4');                                 % teti
xmax = fread(fid,1,'real*4');                                 % xmax
zmax = fread(fid,1,'real*4');                                 % zmax
nnx = fread(fid,1,'integer*4');                               % nnx
nnz = fread(fid,1,'integer*4');                               % nnz

vxs = zeros(nnx,nnz,nss);
vys = zeros(nnx,nnz,nss);
vzs = zeros(nnx,nnz,nss);

for is = 1:nss, vxs(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxs 
for is = 1:nss, vys(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vys 
for is = 1:nss, vzs(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vzs
bx = fread(fid,[nnx nnz],'real*4');                           % bx 
by = fread(fid,[nnx nnz],'real*4');                           % by 
bz = fread(fid,[nnx nnz],'real*4');                           % bz 
ex = fread(fid,[nnx nnz],'real*4');                           % ex 
ey = fread(fid,[nnx nnz],'real*4');                           % ey 
ez = fread(fid,[nnx nnz],'real*4');                           % ez 

dns = zeros(nnx,nnz,nss);
for is = 1:nss, dns(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % dns 
xe = fread(fid,nnx,'real*4');                                 % xe 
ze = fread(fid,nnz,'real*4');                                 % ze 
mass = fread(fid,nss,'real*4');                               % mass 
q = fread(fid,nss,'real*4');                                  % q 
time = fread(fid,1,'real*8');                                 % time 
%time = fread(fid,1,'real*4');                                 % time 
wpewce = fread(fid,1,'real*4');                               % wpewce 
dfac = fread(fid,nss,'real*4');                               % dfac
vxx = zeros(nnx,nnz,nss); % vxvx
vyy = zeros(nnx,nnz,nss); % vyvy
vzz = zeros(nnx,nnz,nss); % vzvz
vxy = zeros(nnx,nnz,nss); % vxvy
vxz = zeros(nnx,nnz,nss); % vxvz
vyz = zeros(nnx,nnz,nss); % vyvz

for is = 1:nss, vxx(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxx 
for is = 1:nss, vyy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vyy 
for is = 1:nss, vzz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vzz 
for is = 1:nss, vxy(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxy 
for is = 1:nss, vxz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vxz 
for is = 1:nss, vyz(:,:,is) = fread(fid,[nnx nnz],'real*4'); end  % vyz
remainder = fread(fid);
st = fclose(fid);

%% Output
% Define all data to be output, two cell arrays, one with char variable
% name and one with data

varstrs = {'it','dt','teti','xmax','zmax','nnx','nnz','xe','ze','mass','q','wpewce','dfac','time',...
  'dns','bx','by','bz','ex','ey','ez','vxs','vys','vzs','vxx','vyy','vzz','vxy','vxz','vyz'};
nvars = numel(varstrs);
vars = cell(1,nvars);

for ivar = 1:nvars
  vars{ivar} = eval(varstrs{ivar});
end

varargout{1} = varstrs;
varargout{2} = vars;