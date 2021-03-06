%% Load data
timestep = 001000;
txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons/parts/parts-%05.0f.dat',timestep); % michael's perturbation

%timestep = 055;
%txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons_2/data/fields-%05.0f.dat',timestep); % try-out with larger perturbation


%tic; [time,r,e,b,n1,ne,ve,vi,je,ji,pe,pi,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q] = read_fields_ieie(txtfile); toc
tic; [x,z,E,B,...
  ni1,ne1,ni2,ne2,...
  vi1,ve1,vi2,ve2,...
  ji1,je1,ji2,je2,...
  pi1,pe1,pi2,pe2,...
  ti1,te1,ti2,te2,...
  dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile); toc