      program read_fields
c
c   read a file, from which the thing can restart
c
      implicit none
c
      include "globalshf.h"
c
      character*5 cdump
      character*5 cdump1
      character*80 ofile
      integer nnx, nnz
c
c  read fields from disk
c
c
c
      if (lead) write(*,*)' reading fields'
      itime=int(time)
      write(cdump,9000)itime
 9000 format(i5.5)
      ofile='data/fields-03000.dat'
      if (lead) write(*,*)ofile
      open(20,file=ofile,form='unformatted')
      read(20)it,dt,teti,xmax,zmax,nnx,nnz,vxs,vys,vzs,
     &         bx,by,bz,ex,ey,ez,dns,xe,ze,mass,q,time,wpewce,dfac,
     &         pxx,pyy,pzz,pxy,pxz,pyz
      close(20)
      if (lead) write(*,*)' fields in'
      print *,it,dt,teti,xmax,zmax,nnx,nnz
c
c
c
c  read particles from disk
c
c      if (lead) write(*,*)' reading particles'
c      if (lead) write(*,*)nsp(1),nsp(2)
c      if (lead) write(*,*)nsact(1),nsact(2)
c      write(cdump,9001)itime
c      write(cdump1,9002)myid
c 9001 format(i5.5)
c 9002 format(i5.5)
c      ofile='parts/parts-'//cdump//'-p'//cdump1//'.dat'
c      open(30,file=ofile,form='unformatted')
c      read(30)it,dt,teti,xmax,zmax,nnx,nnz,x,z,vx,vy,vz,nsp,nsact
c      if (lead) write(*,*)' particles in'
c      close(30)
c      if (lead) write(*,*)nsp(1),nsp(2)
c      if (lead) write(*,*)nsact(1),nsact(2)
c      stop
c
c  assign some stuff that is otherwise left out:
c
c      dx=xe(3)-xe(2)
c      dz=ze(3)-ze(2)
c      if (lead) write(*,*)dx,dz
c
c   the rest is missing, but should not be needed
c     +         bxa,bya,bza,phi,
c     +         tend,
c     +         frac,drive,c,
c     +         emf,l
c
      return
      end
c