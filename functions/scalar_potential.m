function phi = scalar_potential(x,z,ex,ez,ni,ne,varagin)
% Calculates the xz-plane magnetic vector potential.
%   phi = scalar_potential(x,z,ex,ez,ni,ne)
%   grad(phi) = E

if nargin == 1 && strcmp(varargin{1},'plot')
  doPlot = 1;
else
  doPlot = 0;
end
  
% Grid
dx = x(2)-x(1);
dz = z(2)-z(1);
nz = numel(z);
nx = numel(x);
nnx = nx;
nnz = nz;

% Dont put zero right at the edge 
% ixm = 10;
% phi = zeros(nx,nz);
% % Advance up
% phi = -cumsum(Ez,2)*dz;
% phi = phi-cumsum(Ex,1)*dx;
% % Advance to the right
% phi((ixm+1):end,:) = repmat(phi(ixm,:),nx-ixm,1) + -cumsum(Ex((ixm+1):end,:),1)*dx;
% % Advance to the left
% phi((ixm-1):-1:1,:) = repmat(phi(ixm,:),ixm-1,1) + -cumsum(Ex((ixm-1):-1:1,:),1)*dx;

dx = x(2)-x(1);
dz = z(2)-z(1);
odx2=1./dx^2;
odz2=1./dz^2;
fac=1./(2.*(odx2+odz2));

% dns(:,:,1) = dni_h;
% dns(:,:,2) = dne_h;
% dns(:,:,3) = dni;
% dns(:,:,4) = dne;
% chargedensity(1:nnx,1:nnz) = zeros;
% q = [1; -1; 1; -1];

% nss = 4;
% for k=1:nss
%   for iz=2:nnz
%     for ix=2:nnx
% %!RHS of poisson equation is "charge"
%         chargedensity(ix,iz)=chargedensity(ix,iz)+q(k)*dfac(k)*dns(ix,iz,k);
% 
%     end
%   end
% end

%chargedensity = (dni_h+dni)-(dne_h+dne);

%% No initial guess... 
phi = zeros(nnx+1,nnz+1);
for i=2:nnx+1 %phi+ -> phi_i+1
    phi(i,nnz/2) = ...
        phi(i-1,nnz/2) - ex(i-1,nnz/2)*dx;
end
% %E = -grad(phi) -> -Ex = d/dx (phi) -> -Ex = phi+ - phi- /dx
% % phi+ = -Ex + phi-
% 
% %for k=nnz/2+1:nnz
% %   phi(1:nnx,k)=phi(1:nnx,k-1)-dz.*ez(:,k-1);
% %end
% 
for j=(nnz/2)+1:nnz+1
    phi(1:nnx,j) = ...
        phi(1:nnx,j-1) - ez(:,j-1)*dz;
end

for j=(nnz/2):-1:2 % phi- = phi+ + e
    phi(1:nnx,j-1) = ...
        phi(1:nnx,j) + ez(:,j-1)*dz;
end
rhs = phi;
phi = phi.*0;


phi_0 = phi;

%rhs = zeros(nnx+1,nnz+1);
rhs(1:nnx,1:nnz) = ni-ne;
rhs(nnx-1,:) = rhs(2,:); % top and bottom equal 
%rhs(:,nnz-1) = rhs(:,2); % left and right





for jj=1:100
    hhphi=abs(odx2*(circshift(phi,[1 0])+circshift(phi,[-1 0])) -2.*(odx2+odz2)*phi ...
          +odz2*(circshift(phi,[0 1])+circshift(phi,[0 -1]))+rhs);
    hhsigma=abs(ni-ne);
    
    hh=sum(sum(hhphi(2:nnx,2:nnz)))/sum(sum(hhsigma(2:nnx,2:nnz)));

    if hh < 1e-2 
         break
    end
    
     
    for ii=1:50 %shift each row: circshift(m,1,2) , circshift(m,1,1) = column
    
    phi = fac*(odx2*(circshift(phi,[1 0])+ circshift(phi,[-1 0])) + ...
                    odz2*(circshift(phi,[0 1])+ circshift(phi,[0 -1])) + rhs);
                
    phi(1,:) = phi(nnx,:); % first row equal to second last
    phi(nnx+1,:) = phi(2,:); % last row equal to row 2
    
    phi(:,1) = phi(:,2); % first column equal to second
    phi(:,nnz+1) = phi(:,nnz); % last column equal second last
    
    
    
    end
    if doPlot
      imagesc(x,z,phi)
      hca = gca;
      hca.YDir = 'normal';
    end
     fprintf('iteration %g %g \n', jj, hh);
        
end
    
        
    
    
  
     
     








% 
% phi(:) = NaN;
% phi(1,nnz/2)=0;
% 
% 
% % Top and bottom
% % for j=2:nnx
% %    phi(j,nz/2)=phi(j-1,nnz/2)-dx*ex(j-1,nnz/2);
% % end
% 

% 
% 
% 
% phi_ez = zeros(nnx,nnz);
% 
% phi_ez(1:nnx,1:nnz) = -cumsum(ez,2)*dz;
% phi_ex(:,:) = -cumsum(ex)*dx;
% 
% %phi_c = phi_ez + phi_ex;
% 
% %% Use only ez as inital guess... - no start with 0
% phi_c = phi_ex.*0;
% rhs = chargedensity;
% 
% %Set x- and z-boundary to be equal wrt eachother
% phi_c(1,:) = phi_c(end,:);
% phi_c(:,1) = phi_c(:,end);
% 
% % Save initial guess
% phi_0 = phi_c;
% 
% for jj=1:60
%     for ii=1:50 %shift each row: circshift(m,1,2) , circshift(m,1,1) = column
%         
%   0      phi_c = fac*(odx2*(circshift(phi_c,[1 0])+ circshift(phi_c,[-1 0])) + ...
%                     odz2*(circshift(phi_c,[0 1])+ circshift(phi_c,[0 -1])) - rhs);
%                 
% phi_c(:,1) = phi_c(:,2);
% phi_c(:,end) = phi_c(:,end-1);
%     end
%     
%     hh=max(abs(phi_c-phi_0))/max(abs(phi_c));
%      if hh < 1e-2
%          break
%      end
%      
%      phi_0 = phi_c;
%      fprintf('iteration %g %g \n', jj, ii);
% end
% 
% % To check if derivative of phi2 is equal to ez:
% % figure; imagesc(-diff(phi_ez,1,2)'./dz);
% % or 
% % figure; imagesc(-diff(phi_ex,1,1)'./dx);
% 
% for k=nnz/2:-1:1
%    phi(1:nnx-1,k-1)=phi(1:nnx-1,k)+dz*ez(:,k-1);   
% end
% phi(nnx,:)=phi(1,:);
% out=phi(1:nnx-1,1:nnz-1);
% 
end