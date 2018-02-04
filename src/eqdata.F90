! This file is part of GTC version 3 
! GTC version 3 is released under the 3-Clause BSD license:

! Copyright (c) 2002,2010,2016, GTC Team (team leader: Zhihong Lin, zhihongl@uci.edu)
! All rights reserved.

! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, 
!    this list of conditions and the following disclaimer.

! 2. Redistributions in binary form must reproduce the above copyright notice, 
!    this list of conditions and the following disclaimer in the documentation 
!    and/or other materials provided with the distribution.

! 3. Neither the name of the GTC Team nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without 
!    specific prior written permission.

module spline_normalization
  contains
    subroutine normalizeSpline2d(spline,physicalUnit,fluxUnit)
      implicit none
  
      ! arguments
      real,dimension(:,:,:),intent(inout) :: spline
      real,intent(in) :: physicalUnit,fluxUnit
      ! internal variables
      integer :: i,pow
      real :: a

      spline=spline/physicalUnit

      do i = 1,9
        a = 1.0
        do pow = 1,modulo(i-1,3)
          a = a*fluxUnit
        enddo
        spline(i,:,:)=spline(i,:,:)*a
        spline(i,1,:)=spline(i,1,:)/sqrt(a)
      enddo
    end subroutine normalizeSpline2d

    subroutine normalizeSpline1d(spline,physicalUnit,fluxUnit)
      implicit none
  
      ! arguments
      real,dimension(:,:),intent(inout) :: spline
      real,intent(in) :: physicalUnit,fluxUnit
      ! internal variables
      integer :: i,pow
      real :: a
  
      spline=spline/physicalUnit

      do i = 1,3
        a = 1.0
        do pow = 1,i-1
          a = a*fluxUnit
        enddo
        spline(i,:)=spline(i,:)*a
      enddo
    end subroutine normalizeSpline1d
end module spline_normalization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eqdata
  use global_parameters
  use equilibrium
  use particle_array
  use field_array,only:iflux
  implicit none

  integer i,icount,ierror,isp,lemax,lrmax,ndim,nfp,sgn,lsp_iflux
  real spdum(10),rmaj,itemp0,psi_iflux,te_iflux,ne_iflux,r_0,r_1,r_iflux,dpx,dp2
  real,dimension(:,:,:),allocatable :: bcn,bsn,xcn,xsn,zcn,zsn,fcn,fsn
  real,dimension(:,:,:,:),allocatable :: bsptmp
  integer,dimension(:),allocatable :: ntor
  character(len=99) cdum
  real,dimension(:,:,:),allocatable :: tmpsp
  real,external :: sprgpsi,sppsirg
  namelist /physical_parameters/ r0,b0,etemp0,eden0,ulength,utime,rho0,betae,tauee,tauii

! equilibrium and profile are represented using 2D-spline on (psi, theta) or 1D-spline on (psi)
! poloidal flux=[0,psiw], spdpsi=psiw/(lsp-1); theta=[0,2pi), spdtheta=2pi/(lst-1)
! lsp and lst are over-written by spdata.dat when numereq=1
  lsp=100
  lst=101

  ndim=1
! set spline dimensionality
! read spline parameters from spdata.dat
  if(numereq/=0)then
     if(mype==0)then
        isp=110
        open(isp,file='spdata.dat',status='old')
        read(isp,*) cdum
        if (numereq==1) then ! axisymmetric equilibrium from EFIT
           read(isp,102) lsp,lst,lemax,lrmax
           read(isp,101) psiw,ped
           psiw=abs(psiw)
           ped=abs(ped)
           lst=lst+1 !in ORBIT spdtheta=2pi/lst; needs to raise lst by one after reading spdata.
        elseif (numereq==2) then ! non-axisymmetric equilibrium from VMEC
           read(isp,102) lsp,lst
           read(isp,101) psiw
           read(isp,102) ndim,nfp
           sgn=1
	   if(psiw<0)sgn=-1 ! poloidal flux determines the direction of axes
           psiw=sgn*psiw
           ped=psiw
           lst=lst+1
        endif
        write(gtcout,*)"lsp=",lsp,"lst",lst
     endif
     icount=1
     call MPI_BCAST(lsp,icount,mpi_Integer,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(lst,icount,mpi_Integer,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(ndim,icount,mpi_Integer,0,MPI_COMM_WORLD,ierror)
  endif

! allocate memory for equilibrium spline array
  allocate (stpp(lsp),mipp(lsp),mapp(lsp),nepp(3,lsp),tepp(3,lsp),nipp(3,lsp),tipp(3,lsp),&
       zepp(3,lsp),ropp(3,lsp),erpp(3,lsp),qpsi(3,lsp),gpsi(3,lsp),ppsi(3,lsp),rpsi(3,lsp),torpsi(3,lsp),&
       bsp(spdim,lsp,lst),xsp(spdim,lsp,lst),zsp(spdim,lsp,lst),gsp(spdim,lsp,lst),rd(spdim,lsp,lst),&
       tmpsp(spdim,lsp,lst),nu(spdim,lsp,lst),dl(spdim,lsp,lst),spcos(3,lst),spsin(3,lst),rgpsi(3,lsp),cpsi(3,lsp),&
       psirg(3,lsp),psitor(3,lsp),nfpp(3,lsp),tfpp(3,lsp),nfepp(3,lsp),tfepp(3,lsp))

  if(mype==0)then
     if(numereq==0)then 
       !Construct analytic equilibrium and profile assuming high aspect-ratio, concentric cross-section
        call analyticeq
     elseif(numereq==1)then ! use EFIT & TRANSP data
        !EFIT spdata used in GTC: 2D spline b,x,z,rd; 1D spline qpsi,gpsi,torpsi
        !TRANSP prodata used in GTC: tipp,tepp,nepp,zeff,ropp,erpp
         call spdata(isp)
         call prodata
     endif
     if(numereq/=2)then
       !modification for various fielddir in equilibrium
        if(fielddir==1 .or. fielddir==3)then
         qpsi(:,:)=-qpsi(:,:)
         gpsi(:,:)=-gpsi(:,:)
        endif
        if(fielddir==2 .or. fielddir==3)then
          tmpsp=bsp
          do i=1,lst-1
            bsp(1:3,:,i)=tmpsp(7:9,:,lst-i)*spdtheta*spdtheta+tmpsp(4:6,:,lst-i)*spdtheta+tmpsp(1:3,:,lst-i)
            bsp(4:6,:,i)=-(2.0_wp*tmpsp(7:9,:,lst-i)*spdtheta+tmpsp(4:6,:,lst-i))
            bsp(7:9,:,i)=tmpsp(7:9,:,lst-i)
          enddo
          tmpsp=xsp
          do i=1,lst-1
            xsp(1:3,:,i)=tmpsp(7:9,:,lst-i)*spdtheta*spdtheta+tmpsp(4:6,:,lst-i)*spdtheta+tmpsp(1:3,:,lst-i)
            xsp(4:6,:,i)=-(2.0_wp*tmpsp(7:9,:,lst-i)*spdtheta+tmpsp(4:6,:,lst-i))
            xsp(7:9,:,i)=tmpsp(7:9,:,lst-i)
          enddo
          tmpsp=zsp
          do i=1,lst-1
            zsp(1:3,:,i)=tmpsp(7:9,:,lst-i)*spdtheta*spdtheta+tmpsp(4:6,:,lst-i)*spdtheta+tmpsp(1:3,:,lst-i)
            zsp(4:6,:,i)=-(2.0_wp*tmpsp(7:9,:,lst-i)*spdtheta+tmpsp(4:6,:,lst-i))
            zsp(7:9,:,i)=tmpsp(7:9,:,lst-i)
          enddo
          tmpsp=gsp
          do i=1,lst-1
            gsp(1:3,:,i)=tmpsp(7:9,:,lst-i)*spdtheta*spdtheta+tmpsp(4:6,:,lst-i)*spdtheta+tmpsp(1:3,:,lst-i)
            gsp(4:6,:,i)=-(2.0_wp*tmpsp(7:9,:,lst-i)*spdtheta+tmpsp(4:6,:,lst-i))
            gsp(7:9,:,i)=tmpsp(7:9,:,lst-i)
          enddo
        endif
     endif
  endif

!broadcast radial and poloildal spline fit of B field and (X,Z) coordinates for axisymmetric equil.
  if(numereq/=2)then
     icount=spdim*lsp*lst
     call MPI_BCAST(bsp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(xsp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(zsp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  endif

! read data and constructing splines for non-axisymmetric equilibrium from VMEC
  if(numereq==2)then
    allocate(bsptmp(spdim,lsp,lst,mtoroidal),fsp(spdim,lsp,lst))
    if(mype==0)then
      allocate(bcn(lsp,lst,ndim),bsn(lsp,lst,ndim),xcn(lsp,lst,ndim),&
            xsn(lsp,lst,ndim),zcn(lsp,lst,ndim),zsn(lsp,lst,ndim),ntor(ndim),&
            fcn(lsp,lst,ndim),fsn(lsp,lst,ndim))
      call spdata3d(isp,ndim,ntor,bcn,bsn,xcn,xsn,zcn,zsn,fcn,fsn,sgn) ! read equilibrium data
      call prodata ! read profiles
    endif
    icount=spdim*lsp*lst*mtoroidal
    !construct 3D spline
    if(mype==0) call spline3d(lsp,lst,mtoroidal,spdpsi,spdtheta,ndim,ntor,bcn,bsn,bsptmp,sgn)
    call MPI_BCAST(fielddir,1,mpi_Integer,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(bsptmp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
    ! assign each toroidal_domain_location with its own values of spline coeffs
    bsp(:,:,:)=bsptmp(:,:,:,toroidal_domain_location+1)

    if(mype==0) call spline3d(lsp,lst,mtoroidal,spdpsi,spdtheta,ndim,ntor,xcn,xsn,bsptmp,sgn)
    call MPI_BCAST(bsptmp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
    xsp(:,:,:)=bsptmp(:,:,:,toroidal_domain_location+1)

    if(mype==0) call spline3d(lsp,lst,mtoroidal,spdpsi,spdtheta,ndim,ntor,zcn,zsn,bsptmp,sgn)
    call MPI_BCAST(bsptmp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
    zsp(:,:,:)=bsptmp(:,:,:,toroidal_domain_location+1)

    if(mype==0) call spline3d(lsp,lst,mtoroidal,spdpsi,spdtheta,ndim,ntor,fcn,fsn,bsptmp,sgn)
    call MPI_BCAST(bsptmp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
    fsp(:,:,:)=bsptmp(:,:,:,toroidal_domain_location+1)

    deallocate(bsptmp)
    if(mype==0) deallocate(bcn,bsn,xcn,xsn,zcn,zsn,fcn,fsn,ntor)
  endif ! end of constructing splines for non-axisymmetric equilibrium
! End of equilibrium and profile construction

  if(mype==0)then
! Now, define a radial coordinate for radial grid: d(rpsi)/d(rgpsi)=sqrt(T_i)
! Using rgpsi, grid cell size for gather-scatter operations is propotioal to local ion gyroradius
     rgpsi=rpsi !iload=1: uniform marker temperature & uniform MHD
     if(iload > 1)then
        rgpsi(1,1)=0.0
        do i=2,lsp
           rgpsi(1,i)=rgpsi(1,i-1)+(rpsi(1,i)-rpsi(1,i-1))/sqrt(0.5*(tipp(1,i-1)+tipp(1,i)))
        enddo
     endif
     rgpsi(1,:)=rgpsi(1,:)*rpsi(1,lsp)/rgpsi(1,lsp)
     call construct_spline1(lsp,spdpsi,rgpsi)

! inverse spline fit: psitor,psirg
! spline fit cell size
     spdtor=torpsi(1,lsp)/real(lsp-1)
     spdrg=rgpsi(1,lsp)/real(lsp-1)
     psitor(1,1)=0.0 !first point is magnetic axis
     call invert_spline0(0,lsp,spdpsi,spdtor,torpsi,psitor)
     psirg(1,1)=0.0 !first point is magnetic axis
     call invert_spline1(lsp,spdpsi,spdrg,rgpsi,psirg)
  endif

! broadcast equilibrium and proile data to all MPI processes
  icount=10
  spdum(1)=spdpsi
  spdum(2)=spdtheta
  spdum(3)=spdrg
  spdum(4)=spdtor
  spdum(5)=ped
  spdum(6)=psiw
  spdum(7)=r0
  spdum(8)=b0
  spdum(9)=etemp0
  spdum(10)=eden0
  call MPI_BCAST(spdum,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  spdpsi=spdum(1)
  spdtheta=spdum(2)
  spdrg=spdum(3)
  spdtor=spdum(4)
  ped=spdum(5)
  psiw=spdum(6)
  r0=spdum(7)
  b0=spdum(8)
  etemp0=spdum(9)
  eden0=spdum(10)

  spdpsi_inv=1.0/spdpsi
  spdtheta_inv=1.0/spdtheta
  spdrg_inv=1.0/spdrg
  spdtor_inv=1.0/spdtor

! radial spline fit of plasma profiles
  icount=3*lsp
  call MPI_BCAST(nepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(tepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(nipp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(tipp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  if(feload>0)then
    call MPI_BCAST(nfepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(tfepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  else
    call MPI_BCAST(nfpp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(tfpp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  endif
  call MPI_BCAST(zepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(ropp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(erpp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

! radial spline fit of geometry
  call MPI_BCAST(qpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(gpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(cpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(rpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(rgpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(torpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(psirg,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(psitor,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

! poloidal spline fit of sin and cos functions
  icount=3*lst
  call MPI_BCAST(spcos,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(spsin,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

  write(gtcout,*)"psi0=",psi0,", psi1=",psi1,", ped=",ped,", psiw=",psiw

! Normalization
#ifndef _FRC
! Find temp and dens @ iflux 
  r_0=sprgpsi(psi0*ped)/sprgpsi(ped)
  r_1=sprgpsi(psi1*ped)/sprgpsi(ped)
  r_iflux=r_0+(r_1-r_0)*real(iflux)/real(mpsi)     ! r @ iflux
  psi_iflux=sppsirg(r_iflux*sprgpsi(ped))          !psi @ iflux
  isp=max(1,min(lsp-1,ceiling(psi_iflux*spdpsi_inv)))
  dpx=psi_iflux-spdpsi*real(isp-1)
  dp2=dpx*dpx
  ne_iflux=nepp(1,isp)+nepp(2,isp)*dpx+nepp(3,isp)*dp2
  te_iflux=tepp(1,isp)+tepp(2,isp)*dpx+tepp(3,isp)*dp2
  if(mype==0)then
    write(gtcout,*)'psi_iflux = ',psi_iflux/ped
    write(gtcout,*)'te_iflux = ',te_iflux
    write(gtcout,*)'ne_iflux = ',ne_iflux
  endif
! Boost profiles so iflux values become true local values  
  tipp=tipp/te_iflux
  tepp=tepp/te_iflux
  tfpp=tfpp/te_iflux
  tfepp=tfepp/te_iflux
  nipp=nipp/ne_iflux
  nepp=nepp/ne_iflux
  nfpp=nfpp/ne_iflux
  nfepp=nfepp/ne_iflux
  if(inorm==0)then 
     etemp0=te_iflux*etemp0
     eden0=ne_iflux*eden0
  endif
#endif  
! GTC unit of length and time: major radius R_0=1, proton cyclotron frequency omega_i=1
  ulength=r0 !read in major radius: length unit (unit=cm)
  utime=1.0_wp/(9580.0_wp*b0) ! time unit (unit=second) = inverse proton cyclotron frequency
! rho0=rho_s=c_s/omega_0 in GTC unit, c_s=sqrt(T_e0/m_p), omega_0= eB_0/m_p, on-axis T_e0, B_0
  rho0=102.0_wp*sqrt(etemp0)/(b0*ulength)
  tstep=sqrt(aion)*tstep/rho0
!electron beta
  betae=4.03e-11_wp*eden0*etemp0/(b0*b0)
  eta=9.64e-5_wp*eta*eden0*utime

! temperature now uses GTC unit
  tepp=tepp*rho0*rho0
  tipp=tipp*rho0*rho0
  if(feload>0)then
    tfepp=tfepp*rho0*rho0
  else
    tfpp=tfpp*rho0*rho0
  endif
! toroidal rotation and radial gradient of equilibrium electrostatic potential
! ropp- toroidal angular velocty, erpp=-dphi_0dp
  if(numereq/=0)then ! numerical equilibrium unit conversion from TRANSP TO GTC
     erpp=erpp*9.58e9_wp*utime*utime/ulength
     ropp=ropp*utime !convert angular frequency from rad/s to Omega_p (GTC unit)
     !erpp=ropp !radial electric field from the force ballance (Omega_t=-dPhi/dpsi)
  else
     ropp=ropp*rho0*rho0
     erpp=erpp*rho0
     ropp=ropp+erpp
  endif
  if(irotation==1)then ! no rotation
     ropp=0.0
  elseif(irotation==0)then !no rotation & E_r
     ropp=0.0
     erpp=0.0
  endif

! on-axis collision time defined as H&H Eq. 1.17 & Eq.4.24, and P.37 of NRL
  tauee=24.0-log(sqrt(eden0)/etemp0)
#ifdef _FRC
  tauee=2.09e7*(etemp0)**1.5*sqrt(afaste)/(eden0*tauee*utime)
#else
  tauee=2.09e7*(etemp0)**1.5*sqrt(aelectron)/(eden0*tauee*utime)
#endif
  tauei=tauee/sqrt(2.0)
! assume Carbon impurity
  itemp0=etemp0*tipp(1,1)/tepp(1,1)
  tauii=23.0-log(sqrt(6.0*zepp(1,1)*eden0)/itemp0**1.5)
  tauii=2.09e7*(itemp0)**1.5*sqrt(aion)/(eden0*tauii*utime)

  ! copy data to GPU
#ifndef GPU_UM
  !$acc update device(qpsi,gpsi,rgpsi,ropp,erpp,cpsi)
  !$acc update device(bsp,xsp)
#endif

  if(mype==0)then
     write(gtcout,physical_parameters)
     write(gtcout,*)"nue_eff=",qpsi(1,1)*sqrt(aelectron/aion)/(rho0*tauei/zepp(1,1)),&
          "nui_eff=",qpsi(1,1)/(sqrt(tipp(1,1))*tauii/zepp(1,1))
     call FLUSH(gtcout)
  endif

101 format(1p4e18.10)
102 format(6i4)

end subroutine eqdata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spdata(isp)
  use global_parameters
  use equilibrium
  use spline_normalization
  implicit none

  integer i,krip,nrip,isp,ipsi
  real rmaj,d0,brip,wrip,xrip,fx,bnorm,testpsi

! for Boozer coordinates only
     do i=1,lsp
        read(isp,101)bsp(1,i,:) !B-field
        read(isp,101)bsp(2,i,:)
        read(isp,101)bsp(3,i,:)
        read(isp,101)bsp(4,i,:)
        read(isp,101)bsp(5,i,:)
        read(isp,101)bsp(6,i,:)
        read(isp,101)bsp(7,i,:)
        read(isp,101)bsp(8,i,:)
        read(isp,101)bsp(9,i,:)
        read(isp,101)xsp(1,i,:) !X-coordinate
        read(isp,101)xsp(2,i,:)
        read(isp,101)xsp(3,i,:)
        read(isp,101)xsp(4,i,:)
        read(isp,101)xsp(5,i,:)
        read(isp,101)xsp(6,i,:)
        read(isp,101)xsp(7,i,:)
        read(isp,101)xsp(8,i,:)
        read(isp,101)xsp(9,i,:)
        read(isp,101)zsp(1,i,:) !Z-coordinate
        read(isp,101)zsp(2,i,:)
        read(isp,101)zsp(3,i,:)
        read(isp,101)zsp(4,i,:)
        read(isp,101)zsp(5,i,:)
        read(isp,101)zsp(6,i,:)
        read(isp,101)zsp(7,i,:)
        read(isp,101)zsp(8,i,:)
        read(isp,101)zsp(9,i,:)
        read(isp,101)gsp(1,i,:) !Gicobian
        read(isp,101)gsp(2,i,:)
        read(isp,101)gsp(3,i,:)
        read(isp,101)gsp(4,i,:)
        read(isp,101)gsp(5,i,:)
        read(isp,101)gsp(6,i,:)
        read(isp,101)gsp(7,i,:)
        read(isp,101)gsp(8,i,:)
        read(isp,101)gsp(9,i,:)
        read(isp,101)qpsi(:,i) !q
        read(isp,101)gpsi(:,i) !g-current
        read(isp,101)rd(1:3,i,1) !i-current
        read(isp,101)ppsi(:,i) !pressure
        read(isp,101)rpsi(:,i) !minor radius
        read(isp,101)torpsi(:,i) !toroidal flux
!!!FOR David EP case: change psi_tor from negative to positive
        torpsi(:,i)=abs(torpsi(:,i))

        rd(1,i,:)=rd(1,i,1)
        rd(2,i,:)=rd(2,i,1)
        rd(3,i,:)=rd(3,i,1)
     enddo
     rd(4:9,:,:)=0.0
     nu=0.0
     dl=0.0
     read(isp,102) krip,nrip
     read(isp,101) rmaj,d0,brip
     read(isp,101) wrip,xrip
!        do i=1,lsp
!           read(isp,101)ha(1,i,:),ha(2,i,:),ha(3,i,:),ha(4,i,:),ha(5,i,:),ha(6,i,:),ha(7,i,:),&
!            ha(8,i,:),ha(9,i,:)
!           read(isp,101)hb(1,i,:),hb(2,i,:),hb(3,i,:),hb(4,i,:),hb(5,i,:),hb(6,i,:),hb(7,i,:),&
!            hb(8,i,:),hb(9,i,:)
!        enddo
     write(gtcout,*)"ripple=",krip,nrip,rmaj,d0,brip,wrip,xrip
close(isp)
#ifdef _FRC
    bnorm=bsp(1,lsp,1)
#else
    bnorm=bsp(1,1,1)
#endif
  write(gtcout,*)"X(0,0)=",xsp(1,1,1),"B_norm=",bnorm

! make sure all quantities from spdata are normalized to gtc units
! R_0=B_0=1

! poloidal flux values get scaled by B_0*R_0^2
  fx=bnorm*xsp(1,1,1)*xsp(1,1,1)

  psiw=psiw/fx
  ped=ped/fx

  write(gtcout,*) "more eqdata: ped,psiw,rmaj=",ped,psiw,rmaj
  spdtheta=2.0*pi/real(lst-1)
  spdpsi=psiw/real(lsp-1)

! higher order spline terms require additional normalization, 
! since they will be multiplied by normalized psi values...

! normalization arguments are (spline, physical unit, flux unit)
  call normalizeSpline1d(torpsi,fx,fx) ! toroidal flux
  call normalizeSpline1d(qpsi,1.0,fx) ! safety factor, q
  call normalizeSpline1d(gpsi,xsp(1,1,1)*bnorm,fx) ! current
  call normalizeSpline2d(rd,xsp(1,1,1)*bnorm,fx) ! current 
  call normalizeSpline2d(gsp,xsp(1,1,1)/bnorm,fx) ! metric tensor (not used)
  call normalizeSpline1d(ppsi,1.0,fx) ! pressure (not used)
  call normalizeSpline2d(bsp,bnorm,fx) ! magnetic field
  call normalizeSpline1d(rpsi,1.0,fx) ! minor radius (already normalized)
  call normalizeSpline2d(zsp,xsp(1,1,1),fx) ! Z (up-down position)
  call normalizeSpline2d(xsp,xsp(1,1,1),fx) ! X (major radius)
  ! note: it is important that XSP IS NORMALIZED LAST, since
  ! we use the xsp(1,1,1) value to normalize other quantities

! toroidal current used in GTC
  cpsi=rd(1:3,:,1)

! near r=0, r=sqrt(2*q_0*psi), expand r in term of sqrt(psi) in first spline cell only
! spline1 for y=y1+y2*sqrt(psi)+y3*psi
  call construct_spline1(lsp,spdpsi,rpsi)

! verify that B-field and wall flux are consistent!
! compute psiw numerically from B-field
! avoid problems from r=0
  testpsi = 0.5*bsp(1,2,1)/sqrt(1+(qpsi(1,2)*xsp(1,2,1)/(xsp(1,2,1)-1))**2) & !Bpol at this radius
          *pi*(xsp(1,2,1)**2) ! area of this annulus
! all other radial points (with r>0)
  if(psi0>1.0)then
    ipsi=int(lsp/10)+1
  else
    ipsi=lsp
  endif
  do i=2,ipsi-1
    testpsi = testpsi + &
             0.5*(bsp(1,i,1)/sqrt(1+(qpsi(1,i)*xsp(1,i,1)/(xsp(1,i,1)-1))**2)+ & 
                  bsp(1,i+1,1)/sqrt(1+(qpsi(1,i+1)*xsp(1,i+1,1)/(xsp(1,i+1,1)-1))**2)) & !Bpol at this radius
             *pi*(xsp(1,i+1,1)**2 - xsp(1,i,1)**2) ! area of this annulus
  enddo
  do i=ipsi,lsp-1
    testpsi = testpsi + &
             0.5*(bsp(1,i,lst/2)/sqrt(1+(qpsi(1,i)*xsp(1,i,lst/2)/(xsp(1,i,lst/2)-1))**2)+ &
                  bsp(1,i+1,lst/2)/sqrt(1+(qpsi(1,i+1)*xsp(1,i+1,lst/2)/(xsp(1,i+1,lst/2)-1))**2)) & !Bpol at this radius
             *pi*(xsp(1,i+1,lst/2)**2 - xsp(1,i,lst/2)**2) ! area of this annulus
  enddo

  testpsi = testpsi/pi2
  ! check if numerical and user input differ by more than 5%
  if (abs((psiw-testpsi)/psiw) > 0.25)then 
    write(gtcout,*)'WARNING: equilibrium problem'
    write(gtcout,*)'B-field value and wall poloidal flux are inconsistent.'
    write(gtcout,*)'user input psiw=',psiw,'computed from B-field psiw=',testpsi
    write(gtcout,*)'Check normalization.'
    write(gtcout,*)'Note: for input B-field in Tesla, psiw should be in Weber/radian'
    !stop
  else
    write(gtcout,*)'Wall flux and B-field are consistent.'
    write(gtcout,*)'user input psiw=',psiw,'computed from B-field psiw=',testpsi
  endif

101 format(1p4e18.10)
102 format(6i4)

end subroutine spdata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prodata
  use global_parameters
  use equilibrium
  use particle_array,only: qion,qelectron,qfast,qfaste,fload,feload
  implicit none

! lpp = # of quantities
! npp = # of radial grids of profile data (flag transp=0)
!       OR
!       2 + # of radial grids of TRANSP profile data (flag transp=1)
  integer,parameter :: maxrows=5000,lpp=15
  integer :: npp
  integer :: i,j,ioprof,iphysical,transp,ios
  real px,dpx
  character(len=11) :: ppname(lpp)
  character(len=1) :: garbage
  real,allocatable,dimension(:,:) :: ppdata

! read equilibrium spline parameters
  ioprof=220

! read through file once, to determine number of rows of prodata
  open(ioprof,file='profile.dat',status='old')
  read(ioprof,*,iostat=ios) garbage
    do i = 1,maxrows
      read(ioprof,*,iostat=ios) px ! borrow px when counting rows since its real 
      if (ios /= 0) exit
      if (i == maxrows) then
        write(gtcout,*)'Error: Maximum number of rows exceeded while reading profile.dat...'
        write(gtcout,*)'Exiting program...'
        stop
      endif
    enddo
    npp = i-1 ! Subtract extra iteration for failed attempt

    rewind(ioprof)

    ! number of radial points determined, so allocate arrays now
    write(gtcout,*)'Reading',npp,'radial points from profile.dat...'
    allocate(ppdata(lpp,npp))
  
    ! read TRANSP data for plasma profile
    read(ioprof,*)ppname
    write(gtcout,*)"plasma profile data fields:"
    write(gtcout,*)ppname
    read(ioprof,*)ppdata
  close(ioprof)

  transp=0
  if(transp==1)then
    ! extend TRANSP data to magnetic axis
    ppdata(1:3,1) = 0.0
    ppdata(4:10,1)   = ppdata(4:10,2)    
    ppdata(11,1)  = 0.0
    ppdata(12:15,1) = ppdata(12:15,2)

    ! extend TRANSP data to wall
    ppdata(1:4,npp)=ppdata(1:4,npp-1)+0.5*(ppdata(1:4,npp-1)-ppdata(1:,npp-2))
    ppdata(5,npp)=ppdata(3,npp)+ppdata(4,npp)
    ppdata(6:15,npp)=ppdata(6:15,npp-1)+0.5*(ppdata(6:15,npp-1)-ppdata(6:15,npp-2))
  endif

! Make sure poloidal flux is consistent with psi values in spdata.dat.
! This assumes that the last point in profile.dat and the outer value 
! in spdata.dat correspond to each other. Both spdata.dat and profile.dat 
! should use a poloidal flux definition where psi=0 on magnetic axis,
! but may be normalized differently.
  ppdata(1,:)=ppdata(1,:)*psiw/ppdata(1,npp)

! spline grid data using linear interpolation on poloidal flux
  do i=2,lsp-1
     px=real(i-1)*ppdata(1,npp)/real(lsp-1)
     j=1
     do while (px>ppdata(1,j))
         j=j+1
     enddo

     dpx=(ppdata(1,j)-px)/(ppdata(1,j)-ppdata(1,j-1))
     stpp(i)=(1.0-dpx)*ppdata(2,j)+dpx*ppdata(2,j-1)
     mipp(i)=(1.0-dpx)*ppdata(3,j)+dpx*ppdata(3,j-1)
     mapp(i)=(1.0-dpx)*ppdata(5,j)+dpx*ppdata(5,j-1)
     tepp(1,i)=(1.0-dpx)*ppdata(6,j)+dpx*ppdata(6,j-1)
     nepp(1,i)=(1.0-dpx)*ppdata(7,j)+dpx*ppdata(7,j-1)
     tipp(1,i)=(1.0-dpx)*ppdata(8,j)+dpx*ppdata(8,j-1)
     zepp(1,i)=(1.0-dpx)*ppdata(9,j)+dpx*ppdata(9,j-1)
     ropp(1,i)=(1.0-dpx)*ppdata(10,j)+dpx*ppdata(10,j-1)
     erpp(1,i)=(1.0-dpx)*ppdata(11,j)+dpx*ppdata(11,j-1)
     if(feload>0)then
       nfepp(1,i)=(1.0-dpx)*ppdata(14,j)+dpx*ppdata(14,j-1)
       tfepp(1,i)=(1.0-dpx)*ppdata(15,j)+dpx*ppdata(15,j-1)
     else
       nfpp(1,i)=(1.0-dpx)*ppdata(14,j)+dpx*ppdata(14,j-1)
       tfpp(1,i)=(1.0-dpx)*ppdata(15,j)+dpx*ppdata(15,j-1)
     endif
  enddo

! magnetix axis
  stpp(1)=ppdata(2,1)
  mipp(1)=ppdata(3,1)
  mapp(1)=ppdata(5,1)
  tepp(1,1)=ppdata(6,1)
  nepp(1,1)=ppdata(7,1)
  tipp(1,1)=ppdata(8,1)
  zepp(1,1)=ppdata(9,1)
  ropp(1,1)=ppdata(10,1)
  erpp(1,1)=ppdata(11,1)
  if(feload>0)then
    nfepp(1,1)=ppdata(14,1)
    tfepp(1,1)=ppdata(15,1)
  else
    nfpp(1,1)=ppdata(14,1)
    tfpp(1,1)=ppdata(15,1)
  endif

! wall
  stpp(lsp)=ppdata(2,npp)
  mipp(lsp)=ppdata(3,npp)
  mapp(lsp)=ppdata(5,npp)
  tepp(1,lsp)=ppdata(6,npp)
  nepp(1,lsp)=ppdata(7,npp)
  tipp(1,lsp)=ppdata(8,npp)
  zepp(1,lsp)=ppdata(9,npp)
  ropp(1,lsp)=ppdata(10,npp)
  erpp(1,lsp)=ppdata(11,npp)
  if(feload>0)then
    nfepp(1,lsp)=ppdata(14,npp)
    tfepp(1,lsp)=ppdata(15,npp)  
  else
    nfpp(1,lsp)=ppdata(14,npp)
    tfpp(1,lsp)=ppdata(15,npp)
  endif

! end of NSTX profile data

! temperature normalized by T_e, density by n_e on axis
  tipp=tipp/tepp(1,1)
  if(feload>0)then
    tfepp=tfepp/tepp(1,1)
  else
    tfpp=tfpp/tepp(1,1)
  endif
  tepp=tepp/tepp(1,1)
  if(feload>0)then
    nfepp=nfepp/nepp(1,1)
  else
    nfpp=nfpp/nepp(1,1)
  endif
    nepp=nepp/nepp(1,1)
  
! thermal (main) ion density assuming charge neutrality q_e*n_e+q_i*n_i=0
#ifdef _FRC
  if(feload>0)then
    nfepp=nepp
    tfepp=tepp
  endif
  nipp=nepp*abs(qfaste)/abs(qion)
#else
  nipp=nepp*abs(qelectron)/abs(qion)
  if(feload>0)nipp=nipp-(nfepp*qfaste)/qion
  if(fload>0)nipp=nipp-(nfpp*qfast)/qion
#endif

!! smooth profile (WARNING: smoothing does not preserve gradient on ends of profile)
!  do j=1,2
!     tepp(1,2:lsp-1)=0.5*tepp(1,2:lsp-1)+0.25*(tepp(1,1:lsp-2)+tepp(1,3:lsp))
!     nepp(1,2:lsp-1)=0.5*nepp(1,2:lsp-1)+0.25*(nepp(1,1:lsp-2)+nepp(1,3:lsp))
!     tipp(1,2:lsp-1)=0.5*tipp(1,2:lsp-1)+0.25*(tipp(1,1:lsp-2)+tipp(1,3:lsp))
!     nipp(1,2:lsp-1)=0.5*nipp(1,2:lsp-1)+0.25*(nipp(1,1:lsp-2)+nipp(1,3:lsp))
!     tfpp(1,2:lsp-1)=0.5*tfpp(1,2:lsp-1)+0.25*(tfpp(1,1:lsp-2)+tfpp(1,3:lsp))
!     nfpp(1,2:lsp-1)=0.5*nfpp(1,2:lsp-1)+0.25*(nfpp(1,1:lsp-2)+nfpp(1,3:lsp))
!     zepp(1,2:lsp-1)=0.5*zepp(1,2:lsp-1)+0.25*(zepp(1,1:lsp-2)+zepp(1,3:lsp))
!     ropp(1,2:lsp-1)=0.5*ropp(1,2:lsp-1)+0.25*(ropp(1,1:lsp-2)+ropp(1,3:lsp))
!     erpp(1,2:lsp-1)=0.5*erpp(1,2:lsp-1)+0.25*(erpp(1,1:lsp-2)+erpp(1,3:lsp))
!  enddo

!convert -dphi/dr to -dphi/dpsi
  erpp(1,lsp)=0.0
  do i=1,lsp
     erpp(1,i)=erpp(1,i)*rpsi(2,i)
  enddo

! spline fit plasma profile data
  call construct_spline0(0,lsp,spdpsi,tepp)
  call construct_spline0(0,lsp,spdpsi,nepp)
  call construct_spline0(0,lsp,spdpsi,tipp)
  call construct_spline0(0,lsp,spdpsi,nipp)
  if(feload>0)then
    call construct_spline0(0,lsp,spdpsi,tfepp)
    call construct_spline0(0,lsp,spdpsi,nfepp)
  else
    call construct_spline0(0,lsp,spdpsi,tfpp)
    call construct_spline0(0,lsp,spdpsi,nfpp)
  endif
  call construct_spline0(0,lsp,spdpsi,zepp)
  call construct_spline0(0,lsp,spdpsi,ropp)
  call construct_spline0(0,lsp,spdpsi,erpp)

write(gtcout,*)"On-axis electron density=",nepp(1,1),"temperature=",tepp(1,1)
write(gtcout,*)"On-axis ion density=",nipp(1,1),"temperature=",tipp(1,1)
write(gtcout,*)"On-axis felectron density=",nfepp(1,1),"temperature=",tfepp(1,1)
write(gtcout,*)"On-axis fion density=",nfpp(1,1),"temperature=",tfpp(1,1)

! overwrite gtc.in data r0, te0, den0 by using prodata.dat
  iphysical=0
  if(iphysical>0)then
     r0=ppdata(4,1)
     etemp0=ppdata(6,1)
     eden0=1.0e19*ppdata(6,1)
  endif

101 format(1p4e18.10)
102 format(6i4)

end subroutine prodata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Construct analytic equilibrium and profile assuming high aspect-ratio, concentric cross-section
subroutine analyticeq
  use global_parameters
  use equilibrium
  use particle_array,only: qion,qelectron,fload,qfast,feload,qfaste,iload
  implicit none

  logical file_exist
  integer i,j,k,ii,jj
  real pdum,r2(3,lsp),sin2(3,lst),cos2(3,lst),sincos(3,lst),c1,c2,c3,c4,c5,c6,&
       q(3),ti(3),tf(3),tfe(3),te(3),ne(3),ni(3),nf(3),nfe(3),ze(3),&
       itemp0,ftemp0,fden0,fetemp0,feden0,er(3),kti,kni,rfpdelta,rfpalpha,rfpbt,rfpmu,rfpa,rfpb(3,lsp),psihat,&
       rdum,edum,fdum,hdum,r2q2(3,lsp),corr(9,lsp,lst),&
       psiw_analytic,ped_analytic,q_analytic(3),&
       ze_analytic(3),er_analytic(3),itemp0_analytic,ftemp0_analytic,&
       fden0_analytic,ne_analytic(3),te_analytic(3),ti_analytic(3),&
       tf_analytic(3),nf_analytic(3),fetemp0_analytic,feden0_analytic,&
       tfe_analytic(3),nfe_analytic(3),&
       q1,q2,q3

  namelist /equilibrium_parameters/ psiw_analytic,ped_analytic,q_analytic,&
       ze_analytic,er_analytic,itemp0_analytic,ftemp0_analytic,&
       fden0_analytic,ne_analytic,te_analytic,ti_analytic,&
       tf_analytic,nf_analytic,fetemp0_analytic,feden0_analytic,&
       tfe_analytic,nfe_analytic

! Test if the input file gtc.in exists
  inquire(file='gtc.in',exist=file_exist)
  if (file_exist) then
      open(55,file='gtc.in',status='old')
      read(55,nml=equilibrium_parameters)
      close(55)
  else
      write(*,*)'Cannot find file gtc.in !!!'
      stop
  endif
  psiw=psiw_analytic
  ped=ped_analytic
  q(1:3)=q_analytic(1:3)
  ze(1:3)=ze_analytic(1:3)
  er(1:3)=er_analytic(1:3)
  itemp0=itemp0_analytic
  ftemp0=ftemp0_analytic
  fden0=fden0_analytic
  fetemp0=fetemp0_analytic
  feden0=feden0_analytic
  ne(1:3)=ne_analytic(1:3)
  te(1:3)=te_analytic(1:3)
  ti(1:3)=ti_analytic(1:3)
  tf(1:3)=tf_analytic(1:3)
  nf(1:3)=nf_analytic(1:3)
  tfe(1:3)=tfe_analytic(1:3)
  nfe(1:3)=nfe_analytic(1:3)

  if(fload==0)fden0=1.0e-6
  if(feload==0)feden0=1.0e-6

! 1D radial spline for qpsi,gpsi,torpsi,tipp,tepp,nepp,zeff,ropp,erpp
! g-current is 1 to order of epsilon
  gpsi(1,:)=1.0

! c-current is 0 to order of epsilon
  cpsi(1,:)=0.0

  spdpsi=psiw/real(lsp-1)
  spdtheta=2.0*pi/real(lst-1)
  do i=1,lsp
     pdum=spdpsi*real(i-1)/psiw
     torpsi(1,i)=q(1)*pdum*psiw+0.5*q(2)*pdum*pdum*psiw+q(3)*pdum*pdum*pdum*psiw/3.0
     rpsi(1,i)=sqrt(2.0*torpsi(1,i)) ! concentric circular cross-section
     qpsi(1,i)=q(1)+q(2)*pdum+q(3)*pdum*pdum
     if(fieldmodel==1 .or. eqcurrent==1) cpsi(1,i)=rpsi(1,i)*rpsi(1,i)/qpsi(1,i)
     zepp(1,i)=ze(1)+ze(2)*pdum+ze(3)*pdum*pdum
     erpp(1,i)=er(1)+er(2)*pdum+er(3)*pdum*pdum
     nepp(1,i)=1.0+ne(1)*(tanh((ne(2)-pdum)/ne(3))-1.0)
     tepp(1,i)=1.0+te(1)*(tanh((te(2)-pdum)/te(3))-1.0)
     tipp(1,i)=1.0+ti(1)*(tanh((ti(2)-pdum)/ti(3))-1.0)
     if(feload>0)then
       tfepp(1,i)=1.0+tfe(1)*(tanh((tfe(2)-pdum)/tfe(3))-1.0)
       nfepp(1,i)=1.0+nfe(1)*(tanh((nfe(2)-pdum)/nfe(3))-1.0)
     else
       tfpp(1,i)=1.0+tf(1)*(tanh((tf(2)-pdum)/tf(3))-1.0)
       nfpp(1,i)=1.0+nf(1)*(tanh((nf(2)-pdum)/nf(3))-1.0)
     endif
  enddo

! normalization for on-axis quantities
  nepp=nepp/nepp(1,1)
  if(feload>0)then
    nfepp=nfepp/nfepp(1,1)
  else
    nfpp=nfpp/nfpp(1,1)
  endif
  tepp=tepp/tepp(1,1)
  tipp=itemp0*tipp/tipp(1,1)
  if(feload>0)then
    tfepp=fetemp0*tfepp/tfepp(1,1)
    nfepp=feden0*nfepp/nfepp(1,1)
  else
    tfpp=ftemp0*tfpp/tfpp(1,1)
    nfpp=fden0*nfpp/nfpp(1,1)
  endif
! main ion density assuming charge neutrality q_e*n_e+q_i*n_i+q_f*n_f(fe)=0
  if(feload>0)then
    nipp=(nepp*abs(qelectron)-nfepp*qfaste)/qion
  elseif(feload==0 .and. fload>0)then
    nipp=(nepp*abs(qelectron)-nfpp*qfast)/qion
  else
    nipp=nepp*abs(qelectron)/qion
  endif

! 1D spline fit: "spline0(1,..." for y=y1+y2*psi+y3*psi*psi
  call construct_spline0(0,lsp,spdpsi,qpsi)
  call construct_spline0(0,lsp,spdpsi,gpsi)
  call construct_spline0(0,lsp,spdpsi,cpsi)
  call construct_spline0(0,lsp,spdpsi,nepp)
  call construct_spline0(0,lsp,spdpsi,tepp)
  call construct_spline0(0,lsp,spdpsi,tipp)
  call construct_spline0(0,lsp,spdpsi,nipp)
  if(feload>0)then
    call construct_spline0(0,lsp,spdpsi,tfepp)
    call construct_spline0(0,lsp,spdpsi,nfepp)
  else
    call construct_spline0(0,lsp,spdpsi,tfpp)
    call construct_spline0(0,lsp,spdpsi,nfpp)
  endif
  call construct_spline0(0,lsp,spdpsi,torpsi)
  call construct_spline0(0,lsp,spdpsi,erpp)

! near r=0, r=sqrt(2*q_0*psi), expand r in term of sqrt(psi) in first spline cell only
! spline1 for y=y1+y2*sqrt(psi)+y3*psi
  call construct_spline1(lsp,spdpsi,rpsi)

! spline fit for cos and sin functions in domain [0,2pi]: spcos,spsin
  do i=1,lst
     pdum=spdtheta*real(i-1)
     spcos(1,i)=cos(pdum)
     spsin(1,i)=sin(pdum)
  enddo
  call construct_spline0(0,lst,spdtheta,spcos)
  call construct_spline0(0,lst,spdtheta,spsin)

! I-current=c-current in boozer coordinates, rd is not used
  rd=0.0

! r^2 and sin^2 and sin*cos
  do i=1,lsp
     c1=rpsi(1,i)
     c2=rpsi(2,i)
     c3=rpsi(3,i)
     r2(1,i)=c1*c1
     r2(2,i)=2.0*c1*c2
     r2(3,i)=2.0*c1*c3+c2*c2
     q1=qpsi(1,i)
     q2=qpsi(2,i)
     q3=qpsi(3,i)
     r2q2(1,i)=c1**2/q1**2
     r2q2(2,i)=2/q1**2-2*c1**2*q2/q1**3
     r2q2(3,i)=-4/q1**3-4*q2/q1**3+6*c1**2*q2**2/q1**4&
                 -2*c1**2*q3/q1**3
  enddo
  do j=1,lst
     c1=spsin(1,j)
     c2=spsin(2,j)
     c3=spsin(3,j)
     c4=spcos(1,j)
     c5=spcos(2,j)
     c6=spcos(3,j)
     sin2(1,j)=c1*c1
     sin2(2,j)=2.0*c1*c2
     sin2(3,j)=2.0*c1*c3+c2*c2
     sincos(1,j)=c1*c4
     sincos(2,j)=c2*c4+c1*c5
     sincos(3,j)=c1*c6+c2*c5+c3*c4
     cos2(1,j)=c4*c4
     cos2(2,j)=2.0*c4*c5
     cos2(3,j)=2.0*c4*c6+c5*c5
  enddo

! 2D spline functions: b,x,z
! rcos(theta_geo)
  do i=1,lsp
     do j=1,lst
        xsp(1,i,j)=rpsi(1,i)*spcos(1,j)
        xsp(2,i,j)=rpsi(2,i)*spcos(1,j)
        xsp(3,i,j)=rpsi(3,i)*spcos(1,j)
        xsp(4,i,j)=rpsi(1,i)*spcos(2,j)
        xsp(5,i,j)=rpsi(2,i)*spcos(2,j)
        xsp(6,i,j)=rpsi(3,i)*spcos(2,j)
        xsp(7,i,j)=rpsi(1,i)*spcos(3,j)
        xsp(8,i,j)=rpsi(2,i)*spcos(3,j)
        xsp(9,i,j)=rpsi(3,i)*spcos(3,j)
     enddo
  enddo
  if(fieldmodel==1)then
    do i=1,lsp
       do j=1,lst
          xsp(1,i,j)=xsp(1,i,j)-r2(1,i)*sin2(1,j)
          xsp(2,i,j)=xsp(2,i,j)-r2(2,i)*sin2(1,j)
          xsp(3,i,j)=xsp(3,i,j)-r2(3,i)*sin2(1,j)
          xsp(4,i,j)=xsp(4,i,j)-r2(1,i)*sin2(2,j)
          xsp(5,i,j)=xsp(5,i,j)-r2(2,i)*sin2(2,j)
          xsp(6,i,j)=xsp(6,i,j)-r2(3,i)*sin2(2,j)
          xsp(7,i,j)=xsp(7,i,j)-r2(1,i)*sin2(3,j)
          xsp(8,i,j)=xsp(8,i,j)-r2(2,i)*sin2(3,j)
          xsp(9,i,j)=xsp(9,i,j)-r2(3,i)*sin2(3,j)
       enddo
    enddo
  endif

! z=rsin(theta_geo)
  do i=1,lsp
     do j=1,lst
        zsp(1,i,j)=rpsi(1,i)*spsin(1,j)
        zsp(2,i,j)=rpsi(2,i)*spsin(1,j)
        zsp(3,i,j)=rpsi(3,i)*spsin(1,j)
        zsp(4,i,j)=rpsi(1,i)*spsin(2,j)
        zsp(5,i,j)=rpsi(2,i)*spsin(2,j)
        zsp(6,i,j)=rpsi(3,i)*spsin(2,j)
        zsp(7,i,j)=rpsi(1,i)*spsin(3,j)
        zsp(8,i,j)=rpsi(2,i)*spsin(3,j)
        zsp(9,i,j)=rpsi(3,i)*spsin(3,j)
     enddo
  enddo
  if(fieldmodel==1)then
    do i=1,lsp
       do j=1,lst
          zsp(1,i,j)=zsp(1,i,j)+r2(1,i)*sincos(1,j)
          zsp(2,i,j)=zsp(2,i,j)+r2(2,i)*sincos(1,j)
          zsp(3,i,j)=zsp(3,i,j)+r2(3,i)*sincos(1,j)
          zsp(4,i,j)=zsp(4,i,j)+r2(1,i)*sincos(2,j)
          zsp(5,i,j)=zsp(5,i,j)+r2(2,i)*sincos(2,j)
          zsp(6,i,j)=zsp(6,i,j)+r2(3,i)*sincos(2,j)
          zsp(7,i,j)=zsp(7,i,j)+r2(1,i)*sincos(3,j)
          zsp(8,i,j)=zsp(8,i,j)+r2(2,i)*sincos(3,j)
          zsp(9,i,j)=zsp(9,i,j)+r2(3,i)*sincos(3,j)
       enddo
    enddo
  endif

#ifdef _CYLINDER
   ! B=B0 (Cylinder has uniform magnetic field)
     do i=1,lsp
        do j=1,lst
           bsp(1,i,j)=1.0      
           do k=2,9            
              bsp(k,i,j)=0.0
           enddo
        enddo
     enddo
#elif _TOKAMAK
   ! B=1-rcos(theta_geo)+(rcos(theta_geo))**2 to order of epsilon^2
     do i=1,lsp
        do j=1,lst
           bsp(1,i,j)=1.0-xsp(1,i,j)+r2(1,i)*cos2(1,j)
           bsp(2,i,j)=0.0-xsp(2,i,j)+r2(2,i)*cos2(1,j)
           bsp(3,i,j)=0.0-xsp(3,i,j)+r2(3,i)*cos2(1,j)
           bsp(4,i,j)=0.0-xsp(4,i,j)+r2(1,i)*cos2(2,j)
           bsp(5,i,j)=0.0-xsp(5,i,j)+r2(2,i)*cos2(2,j)
           bsp(6,i,j)=0.0-xsp(6,i,j)+r2(3,i)*cos2(2,j)
           bsp(7,i,j)=0.0-xsp(7,i,j)+r2(1,i)*cos2(3,j)
           bsp(8,i,j)=0.0-xsp(8,i,j)+r2(2,i)*cos2(3,j)
           bsp(9,i,j)=0.0-xsp(9,i,j)+r2(3,i)*cos2(3,j)
        enddo
     enddo
     if(fieldmodel==0 .and. iload==9)then
       do i=1,lsp
         do j=1,lst
           corr(1,i,j)=0.5*r2q2(1,i)-0.5*r2(1,i)*spcos(1,j)
           corr(2,i,j)=0.5*r2q2(2,i)-0.5*r2(2,i)*spcos(1,j)
           corr(3,i,j)=0.5*r2q2(3,i)-0.5*r2(3,i)*spcos(1,j)
           corr(4,i,j)=-0.5*r2(1,i)*spcos(2,j)
           corr(5,i,j)=-0.5*r2(2,i)*spcos(2,j)
           corr(6,i,j)=-0.5*r2(3,i)*spcos(2,j)
           corr(7,i,j)=-0.5*r2(1,i)*spcos(3,j)
           corr(8,i,j)=-0.5*r2(2,i)*spcos(3,j)
           corr(9,i,j)=-0.5*r2(3,i)*spcos(3,j)
         enddo
       enddo
       bsp(1:9,:,:)=bsp(1:9,:,:)+corr(1:9,:,:)
     endif
   ! x=1+rcos(theta_geo)
     xsp(1,:,:)=1.0+xsp(1,:,:)
#endif

  ! calculation of the toroidal rotation
  do i=1,lsp
     kni= nipp(2,i)/nipp(1,i)
     kti= tipp(2,i)/tipp(1,i)
     ropp(1,i)=-(kni + (1.0-1.17)*kti)*tipp(1,i)/qion
  enddo
  call construct_spline0(0,lsp,spdpsi,ropp)

end subroutine analyticeq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spdata3d(isp,ndim,ntor,bcn,bsn,xcn,xsn,zcn,zsn,fcn,fsn,sgn)
  use global_parameters
  use equilibrium
  implicit none

  integer i,n,isp

  integer,intent(in) :: ndim,sgn
  integer,intent(out) :: ntor(ndim)
  real(wp),intent(out) :: bcn(lsp,lst,ndim),bsn(lsp,lst,ndim),xcn(lsp,lst,ndim),&
                          xsn(lsp,lst,ndim),zcn(lsp,lst,ndim),zsn(lsp,lst,ndim),&
                          fcn(lsp,lst,ndim),fsn(lsp,lst,ndim)
  real(wp) :: b00,x00,dum


! read toroidal mode numbers
  read(isp,102) ntor(1:ndim)
! read toroidal harmonics of B,X,Z,dPHI=PHI-phi_boozer
  do n = 1, ndim
     do i = 1, lsp
           read(isp,101) bcn(i,1:lst,n)
           read(isp,101) bsn(i,1:lst,n)
           read(isp,101) xcn(i,1:lst,n)
           read(isp,101) xsn(i,1:lst,n)
           read(isp,101) zcn(i,1:lst,n)
           read(isp,101) zsn(i,1:lst,n)
           read(isp,101) fcn(i,1:lst,n)
           read(isp,101) fsn(i,1:lst,n)
     enddo
  enddo
  do i=1,lsp
    read(isp,101)qpsi(1,i) !q
    read(isp,101)gpsi(1,i) !g-current
    read(isp,101)cpsi(1,i) !i-current (plasma current)
    read(isp,101)rpsi(1,i) !minor radius
    read(isp,101)torpsi(1,i) !toroidal flux
  enddo
  close(isp)

! change sign due to VMEC convention for z->-z
  fcn=-fcn
  fsn=-fsn
  torpsi=-torpsi
  cpsi=-cpsi 
  qpsi=-qpsi
 
  if((sgn>0).and.(torpsi(1,lsp)>0))fielddir=0
  if((sgn>0).and.(torpsi(1,lsp)<0))fielddir=1 ! negative q
  if((sgn<0).and.(torpsi(1,lsp)<0))fielddir=2 ! change theta-zeta dir
  if((sgn<0).and.(torpsi(1,lsp)>0))fielddir=3 ! negative q; change theta-zeta dir

  write(gtcout,*)"In zeta out system: sgn(psip)=",sgn,"torpsi(1,lsp)=",torpsi(1,lsp),&
                 "gpsi(1,1)=",gpsi(1,1),'fielddir=',fielddir

! change signs for flipped coordinate system (psip already changed so no change to q)
  torpsi=sgn*torpsi
  gpsi=sgn*gpsi 
  cpsi=sgn*cpsi


  x00=sum(xcn(1,:,1))/lst
  b00=sum(bcn(1,:,1))/lst

  write(gtcout,*)"X(0,0)=",x00,"B(0,0)=",b00,'ndim=',ndim
  call FLUSH(gtcout)

! make sure all quantities from spdata are normalized to gtc units R_0=B_0=1
! poloidal flux values get scaled by B_0*R_0^2

  psiw=psiw/(b00*x00*x00)
  ped=ped/(b00*x00*x00)

  spdtheta=2.0*pi/real(lst-1)
  spdpsi=ped/real(lsp-1)

! toroidal flux also scaled by B_0*R_0^2
  torpsi=torpsi/(b00*x00*x00)

! current normalized by R_0*B_0
  gpsi=gpsi/(b00*x00)
  cpsi=cpsi/(b00*x00)

! Jacobian normalized by R_0*B_0
!  gsp=gsp*b00/x00

! normalize B-field by on-axis value
  bcn=bcn/b00
  bsn=bsn/b00
! normalize spatial vars normalized by major radius
  xcn=xcn/x00
  xsn=xsn/x00

  zcn=zcn/x00
  zsn=zsn/x00

  rpsi=rpsi/x00

  rpsi(1,1)=0.0

! construct 1D poloidal splines
  call construct_spline0(0,lsp,spdpsi,qpsi)
  call construct_spline0(0,lsp,spdpsi,gpsi)
  call construct_spline0(0,lsp,spdpsi,cpsi)
  call construct_spline0(0,lsp,spdpsi,torpsi)
  call construct_spline1(lsp,spdpsi,rpsi)
!  call construct_spline0(0,lsp,spdpsi,rpsi)

  rd=0.0
  nu=0.0
  dl=0.0

101 format(1p4e18.10)
102 format(6i4)

end subroutine spdata3d

