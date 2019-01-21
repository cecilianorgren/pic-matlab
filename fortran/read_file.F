      program read_file

integer nmax, u
      parameter (nmax=1000, u=20)
      real x(nmax), y(nmax), z(nmax)

c  Open the data file
      open (u, FILE='/Users/cno062/Data/PIC/points.dat', STATUS='OLD')

c  Read the number of points
      read(u,*) n
      if (n.GT.nmax) then
         write(*,*) 'Error: n = ', n, 'is larger than nmax =', nmax
         goto 9999
      endif

c  Loop over the data points
      do 10 i= 1, n
         read(u,100) x(i), y(i), z(i)
      10 enddo
      100 format (3(F10.4))

c  Close the file
      close (u)

c  Now we should process the data somehow...
c  (missing part)

 9999 stop
      end