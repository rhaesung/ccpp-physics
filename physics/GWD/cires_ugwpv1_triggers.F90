!>\file cires_ugwpv1_triggers.F90
!!
module cires_ugwpv1_triggers

        use mo_ugwpv1_kind,     only: wp

contains

!
!
!> V1: GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>      
      subroutine slat_geos5_tamp_v1(im, tau_amp, xlatdeg, tau_gw)

      implicit none
      integer :: im     
      real(kind=wp)    :: tau_amp, xlatdeg(im), tau_gw(im)
      real(kind=wp)    :: latdeg, flat_gw, tem
      integer :: i
      
!
! if-lat
!
      do i=1, im
        latdeg = abs(xlatdeg(i))    
        if (latdeg < 15.3_wp) then
          tem = (latdeg-3.0_wp) / 8.0_wp
          flat_gw = 0.75_wp * exp(-tem * tem)
          if (flat_gw < 1.2_wp .and. latdeg <= 3.0_wp) flat_gw = 0.75_wp
        elseif (latdeg <  31.0_wp .and. latdeg >=  15.3_wp) then
           flat_gw =  0.10_wp
        elseif (latdeg <  60.0_wp .and. latdeg >=  31.0_wp) then
          tem = (latdeg-60.0_wp) / 23.0_wp
          flat_gw =  0.50_wp * exp(- tem * tem)
        elseif (latdeg >=  60.0_wp) then
          tem = (latdeg-60.0_wp) / 70.0_wp
          flat_gw =  0.50_wp * exp(- tem * tem)
        endif
        tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
      end subroutine slat_geos5_tamp_v1
!      
!> modified for FV3GFS-127L/C96 QBO-experiments
!! GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch)
      subroutine slat_geos5_2020(im, tau_amp, xlatdeg, tau_gw)

      implicit none
      integer :: im     
      real(kind=wp)    :: tau_amp, xlatdeg(im), tau_gw(im)
      real(kind=wp)    :: latdeg, flat_gw, tem
      real(kind=wp), parameter    :: fampqbo = 1.25_wp    ! 1.5     
      real(kind=wp), parameter    :: famp60S = 1.0_wp     ! 1.5
      real(kind=wp), parameter    :: famp60N = 1.0_wp     ! 1.0
      real(kind=wp), parameter    :: famp30  = 0.25_wp    ! 0.4
        
      real(kind=wp), parameter    :: swid15  = 12.5_wp  
      real(kind=wp), parameter    :: swid60S  = 30.0_wp   ! 40
      real(kind=wp), parameter    :: swid60N  = 25.0_wp   ! 30                 
      integer :: i
!    
!
!
      do i=1, im
      
        latdeg = abs(xlatdeg(i))    
        if (latdeg < 15.3_wp) then
          tem = (latdeg-3.0_wp) / swid15
          flat_gw = fampqbo * exp(-tem * tem)
          if (latdeg <= 3.0_wp) flat_gw = fampqbo
        elseif (latdeg <  31.0_wp .and. latdeg >=  15.3_wp) then
           flat_gw =  famp30
        elseif (latdeg <  60.0_wp .and. latdeg >=  31.0_wp) then
          tem = (latdeg-60.0_wp) / 23.0_wp
          flat_gw =   famp60N* exp(- tem * tem)
        elseif (latdeg >=  60.0_wp) then
          tem = (latdeg-60.0_wp) /swid60N
          flat_gw =  famp60N * exp(- tem * tem)
        endif

          if (xlatdeg(i) <= -31.0_wp) then   
!                 
            if (latdeg <  60.0_wp .and. latdeg >=  31.0_wp) then
              tem = (latdeg-60.0_wp) / 23.0_wp
              flat_gw =  famp60S * exp(- tem * tem)           
            endif           
            if (latdeg >=  60.0_wp) then
               tem = (latdeg-60.0_wp) /swid60S
               flat_gw =  famp60S * exp(- tem * tem)
            endif
            
          endif
        tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
      end subroutine slat_geos5_2020   
      
!>  WAM: GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>           
      subroutine slat_geos5(im, xlatdeg, tau_gw)
      
      implicit none
      integer :: im     
      real(kind=wp)  :: xlatdeg(im)          
      real(kind=wp)  :: tau_gw(im)
      real(kind=wp)  :: latdeg
      real(kind=wp), parameter  :: tau_amp = 3.5e-3_wp    ! 3.5 mPa 
      real(kind=wp)             :: trop_gw,  flat_gw 
      integer :: i
!
! if-lat
!
      trop_gw = 0.75_wp
      do i=1, im
      latdeg = xlatdeg(i)    
      if (-15.3_wp < latdeg .and. latdeg < 15.3_wp) then
          flat_gw = trop_gw*exp(-( (abs(latdeg)-3.0_wp)/8.0_wp)**2)
          if (flat_gw < 1.2_wp .and. abs(latdeg) <= 3.0_wp) flat_gw = trop_gw
      else if (latdeg > -31.0_wp .and. latdeg <= -15.3_wp) then
          flat_gw =  0.10_wp
      else if (latdeg <  31.0_wp .and. latdeg >=  15.3_wp) then
           flat_gw =  0.10_wp
      else if (latdeg > -60.0_wp .and. latdeg <= -31.0_wp) then
        flat_gw =  0.50_wp*exp(-((abs(latdeg)-60.0_wp)/23.0_wp)**2)
      else if (latdeg <  60.0_wp .and. latdeg >=  31.0_wp) then
        flat_gw =  0.50_wp*exp(-((abs(latdeg)-60.0_wp)/23.0_wp)**2)
      else if (latdeg <= -60.0_wp) then
         flat_gw =  0.50_wp*exp(-((abs(latdeg)-60.0_wp)/70.0_wp)**2)
      else if (latdeg >=  60.0_wp) then
         flat_gw =  0.50_wp*exp(-((abs(latdeg)-60.0_wp)/70.0_wp)**2)
      end if
      tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
      end subroutine slat_geos5      
      
!>   Spontaneous GW triggers by dynamical inbalances (OKW, fronts/jets, and convection)
!!       not activated due to "limited" set of GFS-physics 
!!       statein-type ( needs horizontal gradients of winds and temperature, humodity)
      subroutine get_spectra_tau_convgw &
      (nw, im, levs, dcheat,  scheat, precip, icld, xlatd, sinlat, coslat,taub, klev, if_src, nf_src)
!
! temporarily can put GEOS-5/MERRA-2 GW-lat dependent function
!      
      integer                   :: nw, im, levs 
      integer,dimension(im,3)   :: icld   
      real(kind=wp), dimension(im, levs) :: dcheat,  scheat   
      real(kind=wp), dimension(im)       :: precip, xlatd, sinlat, coslat      
      real(kind=wp), dimension(im)       :: taub
      integer, dimension(im)    :: klev, if_src
      integer ::  nf_src
!
! locals
      real(kind=wp), parameter           ::  precip_max = 100.0_wp   ! mm/day
      real(kind=wp), parameter           ::  tau_amp    = 3.5e-3_wp  ! 3.5 mPa   
       
      integer  :: i, k, klow, ktop, kmid
      real(kind=wp)     :: dtot, dmax, daver
!      
      nf_src = 0
      if_src(1:im) = 0 
      taub(1:im)   = 0.0_wp
      do i=1, im
        klow = icld(i,1)
        ktop = icld(i,2)
        kmid=  icld(i,3)
        if (klow == -99 .and.  ktop == -99) then 
          cycle
        else 
          klev(i) = ktop
          k = klow
          klev(i) = k
          dmax = abs(dcheat(i,k) + scheat(i,k))
          do k=klow+1, ktop
            dtot =abs(dcheat(i,k) + scheat(i,k))
            if ( dtot >  dmax) then
              klev(i) = k
              dmax =  dtot
            endif  
          enddo
!         
! klev  as  max( dcheat(i,k) + scheat)
!       vertical width of conv-heating
!  
! counts/triiger=1   & taub(i) 
!
          nf_src = nf_src +1
          if_src(i) = 1          
          taub(i) = tau_amp* precip(i)/precip_max*coslat(i)
        endif 

      enddo
!
!  
!      
      call Slat_geos5(im, xlatd, taub)
      nf_src =im
      do i=1, im
        if_src(i) = 1
        klev(i)   = 127-45
      enddo
    
! with info on precip/clouds/dc_heat create Bulk 
!        taub(im), klev(im)      
!   
!      print *, '   get_spectra_tau_convgw '  
      end subroutine get_spectra_tau_convgw

!>
      subroutine get_spectra_tau_nstgw(nw, im, levs, trig_fgf, xlatd, sinlat, coslat, taub, klev, if_src, nf_src)
      integer                      :: nw, im, levs 
      real(kind=wp), dimension(im, levs)    :: trig_fgf
!      real(kind=wp), dimension(im, levs+1)  :: pint      
      real(kind=wp), dimension(im)          :: xlatd, sinlat, coslat        
      real(kind=wp), dimension(im)          :: taub
      integer, dimension(im)          :: klev, if_src
      integer ::  nf_src
! locals
      real(kind=wp), parameter           ::  tlim_fgf = 100.0_wp     ! trig_fgf > tlim_fgf, launch waves should scale-dependent
      real(kind=wp), parameter           ::  tau_amp   = 3.5e-3_wp  ! 3.5 mPa
      real(kind=wp), parameter           ::  pmax = 750.e2_wp,  pmin = 100.e2_wp  
      integer, parameter        ::  klow =127-92,   ktop=127-45
      integer, parameter        ::  kwidth = ktop-klow+1   
      integer  :: i, k,  kex
      real(kind=wp)     :: dtot, dmax, daver      
      real(kind=wp)     :: fnorm, tau_min     
      nf_src = 0      
      if_src(1:im) = 0 
      taub(1:im)   = 0.0_wp
      fnorm = 1.0_wp / real(kwidth, kind=wp)
      tau_min = tau_amp*fnorm
      do i=1, im
!
! only trop-c fjets so find max(trig_fgf) => klev
!                      use abs-values to scale tau_amp
!
  
        k = klow
        klev(i) = k
        dmax = abs(trig_fgf(i,k))
        kex  = 0
        if (dmax >= tlim_fgf) kex = kex+1 
        do k=klow+1, ktop
          dtot = abs(trig_fgf(i,k))
          if (dtot >= tlim_fgf) kex = kex+1 
          if ( dtot >  dmax) then
            klev(i) = k
            dmax    =  dtot
          endif  
        enddo

        if (dmax .ge. tlim_fgf) then
          nf_src = nf_src +1
          if_src(i) = 1          
          taub(i) = tau_min*real(kex, kind=wp)          !* precip(i)/precip_max*coslat(i)
        endif 

      enddo
!     
!      print *, '   get_spectra_tau_nstgw '
      call Slat_geos5(im, xlatd, taub)
      nf_src =im
      do i=1, im
        if_src(i) = 1
        klev(i) = 127-45                  ! FV3-127L
      enddo
!                 
      end subroutine get_spectra_tau_nstgw   
! 
!>
      subroutine get_spectra_tau_okw(nw, im, levs,  trig_okw, xlatd, sinlat, coslat, taub, klev, if_src, nf_src)
      integer                      :: nw, im, levs 
      real(kind=wp), dimension(im, levs)    :: trig_okw
!      real(kind=wp), dimension(im, levs+1)  :: pint    
      real(kind=wp), dimension(im)          :: xlatd, sinlat, coslat     
      real(kind=wp), dimension(im)          :: taub
      integer, dimension(im)       :: klev, if_src
      integer ::  nf_src
! locals
      real(kind=wp), parameter           ::  tlim_okw = 100.0_wp                  ! trig_fgf > tlim_fgf, launch GWs should scale-dependent
      real(kind=wp), parameter           ::  tau_amp   = 35.e-3_wp               ! 35 mPa  
      real(kind=wp), parameter           ::  pmax = 750.e2_wp,  pmin = 100.e2_wp  
      integer, parameter        ::  klow =127-92,   ktop=127-45      ! for FV3-127L
      integer, parameter        ::  kwidth = ktop-klow+1   
      integer  :: i, k,  kex
      real(kind=wp)     :: dtot, dmax, daver      
      real(kind=wp)     :: fnorm, tau_min       
      
      nf_src = 0  
      if_src(1:im) = 0    
      taub(1:im) = 0.0_wp    
       fnorm = 1.0_wp/real(kwidth, kind=wp)
      tau_min = tau_amp*fnorm 
      print *, '   get_spectra_tau_okwgw '   
      do i=1, im
        k = klow
        klev(i) = k          
        dmax = abs(trig_okw(i,k))
        kex = 0
        if (dmax >= tlim_okw) kex = kex+1 
        do k=klow+1, ktop
          dtot = abs(trig_okw(i,k))
          if (dtot >= tlim_okw ) kex = kex+1 
          if ( dtot >  dmax) then
            klev(i) = k
            dmax =  dtot
          endif  
        enddo
!         
        if (dmax >= tlim_okw) then
          nf_src = nf_src + 1
          if_src(i) = 1          
          taub(i) = tau_min*real(kex, kind=wp)          !* precip(i)/precip_max*coslat(i)
        endif 

      enddo 
      print *, '   get_spectra_tau_okwgw '           
      end subroutine get_spectra_tau_okw   

end module cires_ugwpv1_triggers
