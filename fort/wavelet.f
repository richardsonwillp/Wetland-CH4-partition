c-----
      program stats
c-----
c     spl_num: Number of data
c     dt: time interval of data in second
c     mh2o: Molecular weight of H2O
c     ma: Molecular weight of dry air
c     mco2: Molecular weight of CO2
c     u: Stremwise wind velocity
c     w: Vertical wind velocity
c     t: Sonic temperature
c     q: H2O density
c     c: CO2 densit
c     ch4: CH4 density
c     pres: Atmospheric pressure
c     ta: Air temperature
c     ep: water vapor pressure/atmospheric pressure
c     mq: Mean H2O density
c     mta: Mean air temperature
c     mc: Mean CO2 density
c     mpres: Mean atmospheric pressure
c     rhoa: dry air density
c     rho: moist air density
c     w_: Wavelet coefficients
c     te_: Wavelet covariance for each scale
c     cov_: Covariance
c     cov_xxwv: Covariance calculated from wavelet coefficients
c-----
      implicit none
      integer spl_num
      real dt,mh2o,ma,mco2
      include 'param1'
      parameter ( mh2o = 18.02, ma = 28.964, mco2 = 44.01 )
      integer i,j,ni,nj,maxi
      real u(spl_num),w(spl_num),t(spl_num),q(spl_num),c(spl_num)
      real ch4(spl_num),pres(spl_num),ta(spl_num)
      real ep,mq,mta,mc,mch4,mpres,rhoa,rho,mthum,mch4se
      real qwpl(spl_num),cwpl(spl_num),ch4wpl(spl_num)
      real e(spl_num),pe(spl_num),QQ(spl_num),RRR(spl_num),SS(spl_num)
      real kappa(spl_num),ch4se(spl_num)
      double precision w_w(spl_num),w_t(spl_num),w_q(spl_num)
      double precision w_c(spl_num),w_ch4(spl_num)
      double precision te_wt(20),te_wq(20),te_wc(20),te_wch(20)
      double precision cov_wt,cov_wq,cov_wc,cov_wch
      double precision cov_wtwv,cov_wqwv,cov_wcwv,cov_wchwv

      open(11,file='temp.dat')
      open(21,file='wvlet.txt')

      ni = 0
      do i=1,spl_num
         read(11,*,err=900) u(i),w(i),t(i),q(i),c(i),ch4(i),pres(i),
     +      ta(i)
         ni = ni + 1
      enddo
 900  continue
      write(6,*) "No of data",ni
      nj = anint(log(real(ni))/log(2.0))
      ni = 2**nj

      mq = 0.
      mta = 0.
      mc = 0.
      mch4 = 0.
      mpres = 0.
      do i=1,ni
         mq = mq + q(i)
         mta = mta + ta(i)
         mc = mc + c(i)
         mch4 = mch4 + ch4(i)
         mpres = mpres + pres(i)
      enddo
      mq = mq / real(ni)
      mta = mta / real(ni)
      mc = mc / real(ni)
      mch4 = mch4 / real(ni)
      mpres = mpres / real(ni)

c-----Point-by-point conversion of sonic temperature and gas densities
c-----Based on Detto & Katul (2007, Boundary-Layer Meteorol 122, 205-216),
c-----and Detto et al. (2011, Agric For Meteorol 151, 1312-1324)
      ep = (mq/mh2o) * 8.3143 * (mta+273.15) / (mpres*1000.)
      rhoa = (mpres-ep*mpres)*1000. / 8.3143 / (mta+273.15) * ma
      rho = rhoa + mq
      mthum = 0.
      do i=1,ni
         t(i) = (t(i)+273.15)*(1.-0.514*q(i)/rho)
     +      - 273.15
         mthum = mthum + t(i)
      enddo
      mthum = mthum / real(ni)
      mch4se = 0.
      do i=1,ni
         qwpl(i) = (1+ep/(1.-ep)) * ((q(i)-mq)+mq/(mta+273.15)
     +      *(t(i)-mthum))/mh2o
         qwpl(i) = qwpl(i)*mh2o + mq
         cwpl(i) = c(i)/mco2 +
     +      (mc/mco2)/(rhoa/ma*1000.)*(q(i)-mq)/mh2o*1000. +
     +      (1+ep/(1.-ep))*(mc/mco2)/(mta+273.15)*(t(i)-mthum)
         cwpl(i) = cwpl(i) * mco2
         e(i) = qwpl(i)/mh2o * 8.3143 * (ta(i)+273.15)
         pe(i) = pres(i) * ( 1. + 0.46*(e(i)/1000)/pres(i))
         QQ(i) = - 1.3*10.**(-7.) * ta(i) + 3.7*10.**(-5.)
         RRR(i) = 4.0*10.**(-8.) * ta(i)**2. + 1.1*10.**(-5.) *
     +      ta(i) + 2.18*10.**(-3.)
         SS(i) = 2.0*10.**(-6.) * ta(i)**2. + 9.8*10.**(-4.) *
     +      ta(i) + 0.378
         kappa(i) = QQ(i) * pe(i)**2. + RRR(i) * pe(i) + SS(i)
         ch4se(i) = ch4(i) * kappa(i)
         mch4se = mch4se + ch4se(i)
      enddo
      mch4se = mch4se / real(ni)
      do i=1,ni
         ch4wpl(i) = (ch4se(i)-mch4se) +
     +      mch4/(rhoa/ma*1000.)*(q(i)-mq)/mh2o*1000. +
     +      (1+ep/(1.-ep))*mch4/(mta+273.15)*(t(i)-mthum)
         ch4wpl(i) = ch4wpl(i) + mch4se
      enddo

      call covs(w,t,ni,cov_wt)
      call covs(w,qwpl,ni,cov_wq)
      call covs(w,cwpl,ni,cov_wc)
      call covs(w,ch4wpl,ni,cov_wch)
      write(6,'(a25,e12.5)') "Time integration cov: ",cov_wch

c-----Wavelet analysis
      call wavelet(1,dt,ni,nj,w,t,te_wt,cov_wtwv,w_w,w_t)
      call wavelet(1,dt,ni,nj,w,qwpl,te_wq,cov_wqwv,w_w,w_q)
      call wavelet(1,dt,ni,nj,w,cwpl,te_wc,cov_wcwv,w_w,w_c)
      call wavelet(1,dt,ni,nj,w,ch4wpl,te_wch,cov_wchwv,w_w,w_ch4)
      write(6,'(a25,e12.5)') "Scale integration cov: ",cov_wchwv

      write(21,*) "1",u(1),w_w(1),w_t(1),w_q(1),w_c(1),w_ch4(1)
      j = 1
      maxi = 2
      do i=2,ni
         write(21,*) j,u(i),w_w(i),w_t(i),w_q(i),w_c(i),w_ch4(i)
         if ( i .eq. maxi ) then
            j = j * 2
            maxi = maxi * 2
         endif
      enddo

  990 continue
      close(11)
      close(21)
      stop
      end

c-----
      subroutine covs(x,y,ni,cov)
      implicit none
      integer i,ni
      real x(ni),y(ni)
      double precision cov,meanx,meany

      cov = 0.d0
      meanx = 0.d0
      meany = 0.d0
      do i=1,ni
         cov = cov + x(i) * y(i)
         meanx = meanx + x(i)
         meany = meany + y(i)
      enddo
      meanx = meanx / dble(ni)
      meany = meany / dble(ni)
      cov = cov / dble(ni) - meanx * meany

      return
      end

c-----
      subroutine wavelet(mwave,dt,ni,nj,a,b,t,cov,wa,wb)
c--VARIABLE DESCRIPTION----------------------------------------------
c mwave: index for mother wavelet
c dt: sampling interval
c m: scale index
c wc: wavelet coefficient
c t: wavelet covariance at scale m
c--------------------------------------------------------------------

      implicit none
      integer ni,nj,mwave
      integer i,n,j,m
      real dt
      real a(ni),b(ni)
      double precision wa(ni),wb(ni),wc(nj)
      double precision t(nj),cov

      do i=1,ni
         wa(i) = a(i)
         wb(i) = b(i)
      enddo

c-----------------------------------------------------------
c     fast wavelet transform
c-----------------------------------------------------------
      call wvlet(mwave,ni,wa,0)
      call wvlet(mwave,ni,wb,0)
c-----data format in wa(ni) at this point
c      wa(1)=S13(1)
c      wa(2)=W13(1): wavelet coefficient at scale 13
c      wa(3)=W12(1)
c      wa(4)=W12(2)
c      wa(5)=W11(1)
c      wa(6)=W11(2)
c      wa(7)=W11(3)
c      wa(8)=W11(4)
c      wa(9)=W10(1)
c      ..........
c------------------------

      j = 1
      m = nj
      n = 1
 200  if ( n .lt. ni ) then
         wc(m) = 0.0d0
         do i=n,j
            wc(m) = wc(m) + wa(i+1)*wb(i+1)
         enddo
         t(m) = wc(m) / dble(ni)
         n = n * 2
         j = j + n
         m = m - 1
         goto 200
      endif

      cov = 0.0d0
      do i=1,nj
         cov = cov + t(i)
      enddo

      return
      end

c-----
      subroutine wvlet(mwave,N,W,inv)
c------------------------------------------------------------------------
c mwave: mother wavelet ( 1 = "Haar", 2 = "Daubeches 4", 3 = "Daubeches 6",
c  4 = "Daubeches 8", 5 = "Daubeches 10", 6 = "Least asymmetric 8",
c  6 = "Least asymmetric 8", 8 = "Coiflet 6" )
c  7 = "Daubeches 10" )
c N: number of data
c J: number of scale
c W: data
c inv: if inv=1, then perform inverse transformation
c g: scaling
c h: wavelet
c------------------------------------------------------------------------
c This code was written based on the following references.
c  Okawa, Y., 2005, 波形の特徴抽出のための数学的処理, CQ出版社 (in Japanese)
c  Percival, D. B. and Walden, A. T., 2000, Wavelet Methods for Time Series
c          Analysis. Cambridge University Press.
      implicit none
      integer mwave,N,J,i,L,k,M,level,t,u,ll,mm,nn,uu,inv
      double precision W(N),X(N),g(16),h(16),coef

      J = anint(log(real(N))/log(2.0))
      if ( mwave .eq. 1 ) then
c        write(6,*) 'Haar'
        L = 2
        g(1) = 1.0d0 / sqrt(2.0d0)
        g(2) = g(1)
      elseif ( mwave .eq. 2 ) then
c        write(6,*) 'Daubeches 4'
        L = 4
        coef = 0.125d0 * sqrt(2.0d0)
        g(1) = coef * ( 1.0d0 + sqrt(3.0d0) )
        g(2) = coef * ( 3.0d0 + sqrt(3.0d0) )
        g(3) = coef * ( 3.0d0 - sqrt(3.0d0) )
        g(4) = coef * ( 1.0d0 - sqrt(3.0d0) )
      elseif ( mwave .eq. 3 ) then
c        write(6,*) 'Daubeches 6'
        L = 6
        g(1) = 0.3326705529500825d0
        g(2) = 0.8068915093110924d0
        g(3) = 0.4598775021184914d0
        g(4) = -0.1350110200102546d0
        g(5) = -0.0854412738820267d0
        g(6) = 0.0352262918857095d0
      elseif ( mwave .eq. 4 ) then
c        write(6,*) 'Daubeches 8'
        L = 8
        g(1) = 0.2303778133088964d0
        g(2) = 0.7148465705529154d0
        g(3) = 0.6308807679398587d0
        g(4) = -0.0279837694168599d0
        g(5) = -0.1870348117190931d0
        g(6) = 0.0308413818355607d0
        g(7) = 0.0328830116668852d0
        g(8) = -0.0105974017850690d0
      elseif ( mwave .eq. 5 ) then
c        write(6,*) 'Daubeches 10'
        L = 10
        g(1) = 0.1601023979741929d0
        g(2) = 0.6038292697971895d0
        g(3) = 0.7243085284377726d0
        g(4) = 0.1384281459013203d0
        g(5) = -0.2422948870663823d0
        g(6) = -0.0322448695846381d0
        g(7) = 0.0775714938400459d0
        g(8) = -0.0062414902127983d0
        g(9) = -0.0125807519990820d0
        g(10) = 0.0033357252854738d0
      elseif ( mwave .eq. 6 ) then
c        write(6,*) 'Least asymmetric 8'
        L = 8
        g(1) = -0.107148901418d0 / sqrt(2.0d0)
        g(2) = -0.041910965125d0 / sqrt(2.0d0)
        g(3) = 0.703739068656d0 / sqrt(2.0d0)
        g(4) = 1.136658243408d0 / sqrt(2.0d0)
        g(5) = 0.421234534204d0 / sqrt(2.0d0)
        g(6) = -0.140317624179d0 / sqrt(2.0d0)
        g(7) = -0.017824701442d0 / sqrt(2.0d0)
        g(8) = 0.045570345896d0 / sqrt(2.0d0)
      elseif ( mwave .eq. 7 ) then
c        write(6,*) 'Least asymmetric 10'
        L = 10
        g(1) = 0.038654795955d0 / sqrt(2.0d0)
        g(2) = 0.041746864422d0 / sqrt(2.0d0)
        g(3) = -0.055344186117d0 / sqrt(2.0d0)
        g(4) = 0.281990696854d0 / sqrt(2.0d0)
        g(5) = 1.023052966894d0 / sqrt(2.0d0)
        g(6) = 0.89658164838d0 / sqrt(2.0d0)
        g(7) = 0.023478923136d0 / sqrt(2.0d0)
        g(8) = -0.247951362613d0 / sqrt(2.0d0)
        g(9) = -0.029842499869d0 / sqrt(2.0d0)
        g(10) = 0.027632152958d0 / sqrt(2.0d0)
      elseif ( mwave .eq. 8 ) then
c        write(6,*) 'Coiflet 6'
        L = 6
        g(1) = -0.011070271529d0 * sqrt(2.0d0)
        g(2) = -0.051429972847d0 * sqrt(2.0d0)
        g(3) = 0.272140543058d0 * sqrt(2.0d0)
        g(4) = 0.602859456942d0 * sqrt(2.0d0)
        g(5) = 0.238929728471d0 * sqrt(2.0d0)
        g(6) = -0.051429728471d0 * sqrt(2.0d0)
      endif
      do k=1,L
        h(k) = (-1)**(k-1) * g(L-k+1)
      enddo

      if ( inv .ne. 1 ) then
c-----wavelet transform (Okawa, Y., 2005)
        M = N
        do level=1,J
          M = M / 2
          do t=1,M
            u = 2 * t - 1
            X(M+t) = h(1) * W(u)
            X(t) = g(1) * W(u)
            do k=2,L
              u = u + 1
              if ( u .gt. 2*M ) then
                u = 1
              endif
              X(M+t) = X(M+t) + h(k) * W(u)
              X(t) = X(t) + g(k) * W(u)
            enddo
          enddo
          do t=1,N
            W(t) = X(t)
          enddo
        enddo
      elseif ( inv .eq. 1 ) then
c-----inverse transform
        M = 1
        do level=J,1,-1
          ll = -1
          mm = 0
          do t=1,M
            ll = ll + 2
            mm = mm + 2
            u = t
            uu = u + M
            i = 1
            k = 2
            X(ll) = g(i)*W(u) + h(i)*W(uu)
            X(mm) = g(k)*W(u) + h(k)*W(uu)
            if ( L .gt. 2 ) then
              do nn = 2,L/2
                u = u - 1
                if ( u .le. 0 ) then
                  u = M
                endif
                uu = u + M
                i = i + 2
                k = k + 2
                if ( i .ge. L ) then
                  i = 1
                  k = 2
                endif
                X(ll) = X(ll) + g(i)*W(u) + h(i)*W(uu) 
                X(mm) = X(mm) + g(k)*W(u) + h(k)*W(uu) 
              enddo
            endif
          enddo
          M = 2 * M
          do t=1,M
            W(t) = X(t)
          enddo
        enddo
      endif

      return
      end
