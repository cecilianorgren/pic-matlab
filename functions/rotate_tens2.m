function vargout = rotate_tens(varargin)
% new_tens = R*old_tens*R^T

switch nargin 
  case 2 % constant rotation matrix
    old_xx = varargin{1}.xx;
    old_xy = varargin{1}.xy;
    old_xz = varargin{1}.xz;
    old_yy = varargin{1}.yy;
    old_yz = varargin{1}.yz;
    old_zz = varargin{1}.zz;        
    rx = varargin{2}(1,:); 
    ry = varargin{2}(2,:);
    rz = varargin{2}(3,:);
  case 4 % spatially varying rotation matrix
    old_xx = varargin{1}.xx;
    old_xy = varargin{1}.xy;
    old_xz = varargin{1}.xz;
    old_yy = varargin{1}.yy;
    old_yz = varargin{1}.yz;
    old_zz = varargin{1}.zz;        
    rx = varargin{2};
    ry = varargin{3};
    rz = varargin{4};   
  case 6 % I don't understand this??
    old_xx = varargin{1}.xx;
    old_xy = varargin{1}.xy;
    old_xz = varargin{1}.xz;
    old_yy = varargin{1}.yy;
    old_yz = varargin{1}.yz;
    old_zx = varargin{1}.zz;    
    rx = varargin{4};
    ry = varargin{5};
    rz = varargin{6};   
  case 9
    old_xx = varargin{1};
    old_xy = varargin{2};
    old_xz = varargin{3};
    old_yy = varargin{4};
    old_yz = varargin{5};
    old_zx = varargin{6};    
    rx = varargin{7};
    ry = varargin{8};
    rz = varargin{9};    
  otherwise
    error('Input not recognized.')    
end
%      | rx.x rx.y rx.z |
% r =  | ry.x ry.y ry.z |
%      | rz.x rz.y rz.z |
%
%      | rx.x ry.x rz.x |
% rt = | rx.y ry.y rz.y |
%      | rx.z ry.z rz.z |
%                       
% newT = r*T*rt
%
%        | rx.x rx.y rx.z |     | T.xx T.xy T.xz |     | rx.x ry.x rz.x |
%     =  | ry.x ry.y ry.z | dot | T.yx T.yy T.yz | dot | rx.y ry.y rz.y |
%        | rz.x rz.y rz.z |     | T.zx T.zy T.zz |     | rx.z ry.z rz.z |
%
%        | rx.x*T.xx + rx.y*T.yx + rx.z*T.zx |     | rx.x ry.x rz.x |
%     =  | ry.x*T.xy + ry.y*T.yy + ry.z*T.zy |  dot | rx.y ry.y rz.y |
%        | rz.x*T.xz + rz.y*T.yz + rz.z*T.zz |      | rx.z ry.z rz.z |
%
%        | rTx |'    | rx.x ry.x rz.x |
%     =  | rTy | dot | rx.y ry.y rz.y |
%        | rTz |     | rx.z ry.z rz.z |


%rxt.x = rx.x; rxt.y = ry.x; rxt.z = rz.x;
%ryt.x = rx.y; ryt.y = ry.y; ryt.z = rz.x;
datasize = size(old_xx);
T(1,1,:) = reshape(old_xx,1,prod(datasize)); 
T(2,2,:) = reshape(old_yy,1,prod(datasize)); 
T(3,3,:) = reshape(old_zz,1,prod(datasize)); 
T(1,2,:) = reshape(old_xy,1,prod(datasize)); 
T(2,1,:) = reshape(old_xy,1,prod(datasize)); 
T(1,3,:) = reshape(old_xz,1,prod(datasize)); 
T(3,1,:) = reshape(old_xz,1,prod(datasize)); 
T(2,3,:) = reshape(old_yz,1,prod(datasize)); 
T(3,2,:) = reshape(old_yz,1,prod(datasize)); 

R(1,1,:) = reshape(rx.x,1,prod(datasize));
R(1,2,:) = reshape(rx.y,1,prod(datasize));
R(1,3,:) = reshape(rx.z,1,prod(datasize));
R(2,1,:) = reshape(ry.x,1,prod(datasize));
R(2,2,:) = reshape(ry.y,1,prod(datasize));
R(2,3,:) = reshape(ry.z,1,prod(datasize));
R(3,1,:) = reshape(rz.x,1,prod(datasize));
R(3,2,:) = reshape(rz.y,1,prod(datasize));
R(3,3,:) = reshape(rz.z,1,prod(datasize));


newT = T*0;
% sum over i j
% nT_mn = R_mi * R_nj * T_ij
for mm = 1:3
  for nn = mm:3
    for ii = 1:3
      for jj = 1:3
        newT(mm,nn,:,:) = newT(mm,nn,:,:) + R(mm,ii,:,:).*R(nn,jj,:,:).*T(ii,jj,:,:);
      end
    end
  end
end

new_xx = squeeze(newT(1,1,:,:));
new_xy = squeeze(newT(1,2,:,:));
new_xz = squeeze(newT(1,3,:,:));
new_yy = squeeze(newT(2,2,:,:));
new_yz = squeeze(newT(2,3,:,:));
new_zz = squeeze(newT(3,3,:,:));

% rTx = old_xx.*rx.x + old_xy.*rx.y + old_xz.*rx.z;
% rTy = old_xy.*rx.x + old_yy.*rx.y + old_yz.*rx.z;
% rTz = old_xz.*rx.x + old_yz.*rx.y + old_zz.*rx.z;
% 
% new_xx = rTx.*rx.x + rTx.*rx.y + rTx.*rx.z;
% 
% new_yy = old_xx.*ry.x + old_xy.*ry.y + old_xz.*ry.z;
% new_xy = old_xy.*ry.x + old_yy.*ry.y + old_yz.*ry.z;
% new_xz = old_xz.*rx.x + old_yz.*rx.y + old_zz.*rx.z;
% 
% 
% new_y = old_x.*ry.x + old_y.*ry.y + old_z.*ry.z;
% new_z = old_x.*rz.x + old_y.*rz.y + old_z.*rz.z;

switch nargout
  case 1
    new_tens.xx = new_xx;
    new_tens.xy = new_xy;
    new_tens.xz = new_xz;
    new_tens.yy = new_yy;
    new_tens.yz = new_yz;
    new_tens.zz = new_zz;
    vargout(1) = new_tens;
  case 6
    vargout(1) = new_xx;
    vargout(2) = new_xy;
    vargout(3) = new_xz;
    vargout(4) = new_yy;
    vargout(5) = new_yz;
    vargout(6) = new_zz;
end
