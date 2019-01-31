function out = tensor_product(ax,ay,az,bx,by,bz)
% TENSOR_PRODUCT Calculates the tensor product AB, where A and B are 
%   vectors with three components.
%   out = TENSOR_PRODUCT(Ax,Ay,Ax,Bx,By,Bz)
prod.xx = ax*bx;
prod.xy = ax*by;
prod.xz = ax*bz;
prod.yy = ay*by;
prod.yz = ay*bz;
prod.zz = az*bz;