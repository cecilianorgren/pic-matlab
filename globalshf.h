c     -*- mode: fortran -*-
      integer nx, nz, np, nss
      parameter (nx=3200, nz=1600, np=150000000, nss=4)
      integer ix,iz,nsp(nss),it,nsact(nss)
c
c  nsp(nss) is the number of particles for each species, to
c  total np
c
      real bx(nx,nz),by(nx,nz),bz(nx,nz),
     .     bxa(nx,nz),bya(nx,nz),bza(nx,nz),
     .     ex(nx,nz),ey(nx,nz),ez(nx,nz),
     .     phi(nx+1,nz+1)
      real charge(nx,nz)
      real x(np),z(np),vx(np),vy(np),vz(np)
      real dns(nx,nz,nss),
     .     vxs(nx,nz,nss),vys(nx,nz,nss),
     .     vzs(nx,nz,nss)
      real xmax,zmax,dt,dx,dz,teti
      real*8 time,tend,toutf,toutp
      real frac(nss),dfac(nss),xmin,zmin
      real zzero
      real xe(nx),ze(nz),mass(nss),q(nss),c,l
      logical drive
      real    wpewce, emf, zBy
      real pxx(nx,nz,nss), pyy(nx,nz,nss), pzz(nx,nz,nss),
     .     pxy(nx,nz,nss), pxz(nx,nz,nss), pyz(nx,nz,nss)
      real qxxx(nx,nz,nss), qxxy(nx,nz,nss), qxxz(nx,nz,nss),
     .     qxyy(nx,nz,nss), qxyz(nx,nz,nss), qxzz(nx,nz,nss),
     .     qyyy(nx,nz,nss), qyyz(nx,nz,nss), qyzz(nx,nz,nss),
     .     qzzz(nx,nz,nss)
      common /heatfl/qxxx,qxxy,qxxz,qxyy,qxyz,qxzz,qyyy,qyyz,qyzz,qzzz
c
c
      integer numprocs,myid,ierr
c	,status(MPI_STATUS_SIZE)
      logical lead, parallel
      real buffers(nx,nz),bufferl(nx,nz,nss)
      integer i, lenl, lens
      integer itime
      common /chat/ numprocs,myid,ierr,parallel,lead,lenl,lens
      common /pressure/  pxx,pyy,pzz,pxy,pxz,pyz
      common /particles/ x,z,vx,vy,vz,nsp,nsact
      common /moments/   dns,vxs,vys,vzs
      common /fields/    bx,by,bz,bxa,bya,bza,ex,ey,ez,phi,charge
      common /params/    time,xmax,zmax,dt,teti,tend,
     .                   frac,dfac,drive,mass,q,c,it,
     .                   zzero
      common /grid/      xe,ze,dx,dz
      common /scale/     wpewce, emf
      common /equil/     l
      common /GuideBC/   zBy
c
c   write out particles over this range only
c
      real xminp, xmaxp, zminp, zmaxp
      common /pwrite/xminp,xmaxp,zminp,zmaxp
c
