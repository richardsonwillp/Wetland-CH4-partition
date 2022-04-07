c-----
      program filter
c-----
      implicit none
      integer ni
      real zm,lb,ub,sd_tch,sd_qch,sd_cch,sd_qt
      include 'param2'
      include 'param3'
      include 'param4'
      include 'param5'
      include 'param6'
      include 'param7'
      integer i,j,k,nj,mj,num,filt_tch,filt_qch,filt_nl,filt_lfc,index
      integer filt_cch
      real date,ws,wd,ts,frq
      real zc,d,sigma_CH4
      real slope_tch,slope_qch,slope_cch,slope_qt,u
      real RMSD_tch,RMSD_qch,RMSD_cch
      real dt
      double precision w,T,q,CO2,CH4,wCH4org,wCH4filt_tch,wCH4filt_qch
      double precision wCH4filt_cch,wCH4filt_lfc,wCH4filt_nl
      character dum
      
      open(11,file='temp2.txt')
      open(12,file='winddata')
      open(13,file='out1.txt')
      open(14,file='out2.txt')
      open(15,file='out3.txt')
      open(16,file='out4.txt')
      open(21,file='filter.txt')
      open(22,file='part.txt')

c-----Load Canopy height, sigma_CH4 data, and DT values      
      open(19,file='zc.txt')
      read(19,*) zc
      open(18,file='sigma_CH4.txt')
      read(18,*) sigma_CH4
      open(17,file='DT.txt')
      read(17,*) dt

c-----multiply sd_xch values by sigma_CH4 to get RMSD defining period's
c-----ebullition threshold
      RMSD_tch = sd_tch*sigma_CH4
      RMSD_qch = sd_qch*sigma_CH4
      RMSD_cch = sd_cch*sigma_CH4

      mj = anint(log(real(ni))/log(2.0))
      d = 0.67 * zc
      read(11,*) dum
      read(11,*) index,u,w,T,q,CO2,CH4
      write(21,*) 'index  w  T  q  CO2  CH4  bubble_index_T-CH4',
     +  '  bubble_index_H2O-CH4 bubble_index_CO2_CH4',
     +  '  lower_frequency_components_index'
      write(21,*) "1",w,T,q,CO2,CH4,"0 0 0 1"
      read(12,*) ws
      read(13,*) slope_qt
      read(14,*) slope_tch
      read(15,*) slope_qch
      read(16,*) slope_cch

      k = 0
      wCH4org = 0.
      wCH4filt_tch = 0.
      wCH4filt_qch = 0.
      wCH4filt_cch = 0.
c-----For lower frequency components
      wCH4filt_lfc = 0.
c-----For non-local components
      wCH4filt_nl = 0.
      nj = 1
      do i=1,mj
         do j=1,nj
            read(11,*) index,u,w,T,q,CO2,CH4
            wCH4org = wCH4org + w*CH4
            ts = real(ni) / real(nj) * dt 
            frq  = (1./ts) * (zm - d) / ws
            if ( frq .ge. lb ) then
               if ( CH4.gt.slope_tch*T+3.*RMSD_tch .or.
     +            CH4.lt.slope_tch*T-3.*RMSD_tch ) then
c-----------------Bubble components
                  filt_tch = 1
               else
                  filt_tch = 0
               endif
               if ( CH4.gt.slope_qch*q+3.*RMSD_qch .or.
     +            CH4.lt.slope_qch*q-3.*RMSD_qch ) then
                  filt_qch = 1
               else
                  filt_qch = 0
               endif
               if ( CH4.gt.slope_cch*CO2+3.*RMSD_cch .or.
     +            CH4.lt.slope_cch*CO2-3.*RMSD_cch ) then
                  filt_cch = 1
               else
                  filt_cch = 0
               endif
               filt_lfc = 0
            else
               filt_tch = 0
               filt_qch = 0
               filt_cch = 0
               filt_lfc = 1
               wCH4filt_lfc = wCH4filt_lfc + w*CH4
            endif
c-----------filtering for non-local processes
c            if ( T.gt.slope_qt*q+3.*sd_qt .or.
c     +         T.lt.slope_qt*q-3.*sd_qt ) then
c               filt_tch = 1
c               filt_qch = 1
c               filt_nl = 1
c               wCH4filt_nl = wCH4filt_nl + w*CH4
c               k = k + 1
c            else
c               filt_nl = 0
c            endif
            wCH4filt_tch = wCH4filt_tch + w*CH4*filt_tch
            wCH4filt_qch = wCH4filt_qch + w*CH4*filt_qch
            wCH4filt_cch = wCH4filt_cch + w*CH4*filt_cch
            write(21,*) nj,w,T,q,CO2,CH4,filt_tch,filt_qch,filt_cch,
     +         filt_lfc
         enddo
         nj = nj * 2
      enddo
      wCH4org = wCH4org / 2.**real(mj)
      wCH4filt_tch = wCH4filt_tch / 2.**real(mj)
      wCH4filt_qch = wCH4filt_qch / 2.**real(mj)
      wCH4filt_cch = wCH4filt_cch / 2.**real(mj)
      wCH4filt_lfc = wCH4filt_lfc / 2.**real(mj)
      wCH4filt_nl = wCH4filt_nl / 2.**real(mj)
      write(22,*) wCH4org,wCH4filt_tch,wCH4filt_qch,wCH4filt_cch,
     +   wCH4filt_lfc
      write(6,'(a31,e13.4)') "Total CH4 flux (mmol m-2 s-1): ",wCH4org
      write(6,'(a46,e13.4)')
     +   "Diffusive CH4 flux (mmol m-2 s-1) based on T: ",
     +   wCH4org-wCH4filt_tch-wCH4filt_lfc
      write(6,'(a47,e13.4)')
     +   "Ebullitive CH4 flux (mmol m-2 s-1) based on T: ",wCH4filt_tch
      write(6,'(a46,e13.4)')
     +   "Diffusive CH4 flux (mmol m-2 s-1) based on q: ",
     +   wCH4org-wCH4filt_qch-wCH4filt_lfc
      write(6,'(a47,e13.4)')
     +   "Ebullitive CH4 flux (mmol m-2 s-1) based on q: ",wCH4filt_qch
      write(6,'(a46,e13.4)')
     +   "Diffusive CH4 flux (mmol m-2 s-1) based on c: ",
     +   wCH4org-wCH4filt_cch-wCH4filt_lfc
      write(6,'(a47,e13.4)')
     +   "Ebullitive CH4 flux (mmol m-2 s-1) based on c: ",wCH4filt_cch
      write(6,'(a40,e13.4)') "Low-frequency component (mmol m-2 s-1): ",
     +   wCH4filt_lfc

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(21)
      close(22)
      stop
      end

