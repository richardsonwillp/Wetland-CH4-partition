c-----
      program ranerror
c-----
      implicit none
      integer i,ni
      include 'param6'
      integer index(ni),filt_tch,filt_qch,filt_lfc,filt_cch
      double precision w(ni),T(ni),q(ni),CH4(ni),CO2(ni)
      double precision CH4dif_tch(ni),CH4bub_tch(ni),CH4dif_qch(ni)
      double precision CH4bub_qch(ni),CH4dif_cch(ni),CH4bub_cch(ni)
      double precision re_wCH4dif_tch,re_wCH4bub_tch,re_wCH4dif_qch
      double precision re_wCH4bub_qch,re_wCH4dif_cch,re_wCH4bub_cch
      character dum

      open(11,file='filter.txt')
      open(21,file='error.txt')

      read(11,*) dum
      do i=1,ni
         read(11,*) index(i),w(i),T(i),q(i),CO2(i),CH4(i),filt_tch,
     +      filt_qch,filt_cch,filt_lfc
         CH4dif_tch(i) = CH4(i)*(1-filt_lfc)*(1-filt_tch)
         CH4bub_tch(i) = CH4(i)*(1-filt_lfc)*filt_tch
         CH4dif_qch(i) = CH4(i)*(1-filt_lfc)*(1-filt_qch)
         CH4bub_qch(i) = CH4(i)*(1-filt_lfc)*filt_qch
         CH4dif_cch(i) = CH4(i)*(1-filt_lfc)*(1-filt_cch)
         CH4bub_cch(i) = CH4(i)*(1-filt_lfc)*filt_cch
         
      enddo

      call wvlet(1,ni,w,1)
      call wvlet(1,ni,CH4dif_tch,1)
      call wvlet(1,ni,CH4bub_tch,1)
      call wvlet(1,ni,CH4dif_qch,1)
      call wvlet(1,ni,CH4bub_qch,1)
      call wvlet(1,ni,CH4dif_cch,1)
      call wvlet(1,ni,CH4bub_cch,1)

      call rferror(ni,w,CH4dif_tch,200,re_wCH4dif_tch)
      call rferror(ni,w,CH4bub_tch,200,re_wCH4bub_tch)
      call rferror(ni,w,CH4dif_qch,200,re_wCH4dif_qch)
      call rferror(ni,w,CH4bub_qch,200,re_wCH4bub_qch)
      call rferror(ni,w,CH4dif_cch,200,re_wCH4dif_cch)
      call rferror(ni,w,CH4bub_cch,200,re_wCH4bub_cch)

      write(21,*) sqrt(re_wCH4dif_tch),sqrt(re_wCH4bub_tch),
     +   sqrt(re_wCH4dif_qch),sqrt(re_wCH4bub_qch),sqrt(re_wCH4dif_cch),
     +   sqrt(re_wCH4bub_cch)

      write(6,'(a55,e13.4)')
     +   'Random diffusive flux error (mmol m-2 s-1) based on T: ',
     +   sqrt(re_wCH4dif_tch)
      write(6,'(a56,e13.4)')
     +   'Random ebullitive flux error (mmol m-2 s-1) based on T: ',
     +   sqrt(re_wCH4bub_tch)
      write(6,'(a55,e13.4)')
     +   'Random diffusive flux error (mmol m-2 s-1) based on q: ',
     +   sqrt(re_wCH4dif_qch)
      write(6,'(a56,e13.4)')
     +   'Random ebullitive flux error (mmol m-2 s-1) based on q: ',
     +   sqrt(re_wCH4bub_qch)
      write(6,'(a55,e13.4)')
     +   'Random diffusive flux error (mmol m-2 s-1) based on c: ',
     +   sqrt(re_wCH4dif_cch)
      write(6,'(a56,e13.4)')
     +   'Random ebullitive flux error (mmol m-2 s-1) based on c: ',
     +   sqrt(re_wCH4bub_cch)

      close(11)
      close(21)
      stop
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

c-----
      subroutine rferror(ni,x,y,m,var)
      implicit none
      integer ni,m
      double precision x(ni),y(ni),var
      integer i,j
      real meanx,meany,gam_xx(m+1),gam_yy(m+1),gam_xy(m+1),gam_yx(m+1)

      meanx = 0.
      meany = 0.
      do i=1,ni
         meanx = meanx + x(i)
         meany = meany + y(i)
      enddo
      meanx = meanx / real(ni)
      meany = meany / real(ni)

      do j=1,m+1
         gam_xx(j) = 0.
         gam_yy(j) = 0.
         gam_xy(j) = 0.
         gam_yx(j) = 0.
         do i=1,ni-j+1
            gam_xx(j) = gam_xx(j) + (x(i)-meanx) * (x(i+j-1)-meanx)
            gam_yy(j) = gam_yy(j) + (y(i)-meany) * (y(i+j-1)-meany)
            gam_xy(j) = gam_xy(j) + (x(i)-meanx) * (y(i+j-1)-meany)
            gam_yx(j) = gam_yx(j) + (y(i)-meany) * (x(i+j-1)-meanx)
         enddo
         gam_xx(j) = gam_xx(j) / real(ni)
         gam_yy(j) = gam_yy(j) / real(ni)
         gam_xy(j) = gam_xy(j) / real(ni)
         gam_yx(j) = gam_yx(j) / real(ni)
      enddo

      var = gam_xx(1) * gam_yy(1) + gam_xy(1) * gam_yx(1)
      do j=2,m+1
         var = var + gam_xx(j)*gam_yy(j)*2. + gam_xy(j)*gam_yx(j)*2.
      enddo
      var = var / real(ni)

      return
      end
