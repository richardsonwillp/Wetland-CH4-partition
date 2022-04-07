c-----
      program pickup
c-----
      implicit none
      integer ni
      real zm,lb,ub
      include 'param2'
      include 'param3'
      include 'param6'
      integer i,j,nj,num,k,mj,index(ni)
      real date,ws,wd,ts,frq,u(ni)
      real zc,d
      real dt
      double precision w(ni),T(ni),q(ni),CH4(ni),CO2(ni)
      character dum

      open(11,file='temp2.txt')
      open(12,file='winddata')
      open(21,file='temp3.txt')
      
      open(19,file='zc.txt')
      read(19,*) zc
      open(17,file='DT.txt')
      read(17,*) dt

      d = 0.67 * zc

      mj = anint(log(real(ni))/log(2.0))
      read(11,*) dum
      ws = 0.
      do i=1,ni
         read(11,*) index(i),u(i),w(i),T(i),q(i),CO2(i),CH4(i)
         ws = ws + u(i)
      enddo
      ws = ws / real (ni)
      write(12,*) ws

      k = 1
      write(21,*) 'index  w  T  q  CO2  CH4'
      nj = 1
      do i=1,mj
         do j=1,nj
            k = k + 1
            ts = real(ni) / real(nj) * dt 
            frq = (1./ts) * (zm - d) / ws
            if ( frq .ge. lb .and. frq .le. ub ) then
               write(21,*) index(k),w(k),T(k),q(k),CO2(k),CH4(k)
            endif
         enddo
         nj = nj * 2
      enddo

      close(11)
      close(12)
      close(21)
      close(19)
      close(17)
      stop
      end

