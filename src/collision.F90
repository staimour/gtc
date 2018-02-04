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

subroutine lorentz(species_name)
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use precision
  implicit none
  
  character(*),intent(in) :: species_name
  
  interface
    subroutine collision_pitch(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                               ncyclep)
      use precision
      !declaration of the dummy arguments
      integer mp,pcoll
      real(wp) qpart,apart,taupp
      real(wp),dimension(0:) :: mesht,meshn
      real(wp),dimension(:,:) :: tppp
      real(wp),dimension(:,:) :: zpart
      integer,optional :: ncyclep
    end subroutine collision_pitch
  end interface

  !call collision_pitch(zion,meshti,meshni,qion,aion,tipp,mi,icoll,tauii)
  if(species_name=="fast-electron")then
    call collision_pitch(zfaste,meshtfe,meshnfe,qfaste,afaste,tfepp,mfe,ecoll,tauei)!,ncyclefe)
  elseif(species_name=="thermal-electron")then  
    call collision_pitch(zelectron,meshte,meshne,qelectron,aelectron,tepp,me,ecoll,tauei)!,ncyclee)  
  endif

end subroutine lorentz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine collision_pitch(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                        ncyclep)
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  implicit none

  !declaration of the dummy arguments
  integer mp,pcoll
  real(wp) qpart,apart,taupp
  real(wp),dimension(0:) :: mesht,meshn
  real(wp),dimension(:,:) :: tppp
  real(wp),dimension(:,:) :: zpart
  integer,optional :: ncyclep
  
  !declaration of the local arguments
  integer m,isp,jst,ii
  real(8) v,dele,zv,ap_inv,ve,nuei,freq,upara,vperp2,pitch,b,randome(mp),rg,dp1,&
       psitmp,dpx,dp2,thetatmp,dtx,dt2,zeff,eden,etemp,ve0,ptmp,tp_inv,delr,delp(mpsi)

  ap_inv=1.0/apart
  tp_inv=1.0/tppp(1,1)
  delr=1.0/deltar
  delp=1.0/deltap
  ve0=rho0*sqrt(2.0*ap_inv)
  
  if(present(ncyclep))then
    nuei=1.88*real(pcoll)*tstep/(taupp*2.0*real(ncyclep))
  else
    nuei=1.88*real(pcoll)*tstep/(taupp*2.0)
  endif
  
  call random_number(randome(1:mp))

!$omp parallel do private(m,psitmp,isp,dpx,dp2,thetatmp,jst,dtx,dt2,rg,ii,dp1,zeff,eden,etemp,b,&
!$omp& upara,vperp2,v,pitch,ve,zv,freq,ptmp)
  do m=1,mp

! radial spline grid
     psitmp=zpart(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
!     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline grid
     thetatmp=zpart(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

! 1D spline in radial
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     zeff=dp1*meshze(ii-1)+(1.0-dp1)*meshze(ii)
     eden=dp1*meshn(ii-1)+(1.0-dp1)*meshn(ii)
     etemp=(dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii))*tp_inv

! 2D spline in (psi, theta) for B-field
     b= bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
       (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
       (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
     
     upara=zpart(4,m)*b*abs(qpart)*ap_inv
     vperp2=zpart(6,m)*zpart(6,m)*2.0*b*ap_inv
     v=sqrt(upara*upara+vperp2)
     pitch=upara/v
     
     ve=ve0*sqrt(etemp)
     zv=max(0.1,min(10.0,v/ve))
     
! nui_ei for electron-ion collision
     freq=nuei*zeff*eden/(zv*zv*zv*etemp**1.5)

! uniform square Montle-Carlo method
     pitch=pitch*(1.0-freq)+(randome(m)-0.5)*sqrt(12.0*freq*max(1.0e-10,1.0-pitch*pitch))
     ptmp=aint(pitch)
     if(abs(ptmp) > 0.5)pitch=sign(1.0,pitch)-(pitch-ptmp)

     upara=v*pitch
     zpart(4,m)=upara*apart/(abs(qpart)*b)
     vperp2=v*v-upara*upara
     zpart(6,m)=sqrt(0.5*apart*(vperp2)/b)
  enddo
  
  return
end subroutine collision_pitch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fokkerplanck(species_name)
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use precision
  implicit none
  character(*),intent(in) :: species_name
  
  interface
    subroutine collision_fokkerplanck(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                                        ncyclep,z_eff)
        use precision
        !declaration of the dummy arguments
        integer mp,pcoll
        real(wp) qpart,apart,taupp
        real(wp),dimension(0:) :: mesht,meshn
        real(wp),dimension(:,:) :: tppp
        real(wp),dimension(:,:) :: zpart
        integer,optional :: ncyclep,z_eff
    end subroutine collision_fokkerplanck
  end interface
  
  if(species_name=='fast-electron')then
#ifdef _FRC
    call collision_fokkerplanck(zfaste,meshtfe,meshnfe,qfaste,afaste,tfepp,mfe,ecoll,tauee)
#else 
    call collision_fokkerplanck(zfaste,meshtfe,meshnfe,qfaste,afaste,tfepp,mfe,ecoll,tauee)!,ncyclep=ncyclefe)
#endif
  elseif(species_name=='thermal_electron')then
    call collision_fokkerplanck(zelectron,meshte,meshne,qelectron,aelectron,tepp,me,ecoll,tauee)!,ncyclep=ncyclee)  
  elseif(species_name=='thermal_ion')then
    call collision_fokkerplanck(zion,meshti,meshni,qion,aion,tipp,mi,icoll,tauii,z_eff=1)
  endif

end subroutine fokkerplanck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine collision_fokkerplanck(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                                  ncyclep,z_eff)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  implicit none

  !declaration of the dummy arguments
  integer mp,pcoll
  real(wp) qpart,apart,taupp
  real(wp),dimension(0:) :: mesht,meshn
  real(wp),dimension(:,:) :: tppp
  real(wp),dimension(:,:) :: zpart
  integer,optional :: ncyclep,z_eff

  integer m,k,icount,ierror,ip,jt,kz,mcell(mp),isp,jst,ii
  real(wp) upara,vperp2,delu,delv2,dneop,dneot,dneoz,psitmp,thetatmp,zetatmp,dpx,dp2,dtx,dt2,b,&
       v,xv,zv,vmin,freq,phix,dphi,f,g,h,sp,sn,dp,dn,dpn,ap_inv,vp,vp0,nupp,pden,ptemp,zeff,&
       random1(mp),random2(mp),bgc(mp),delm(neop*neot*neoz),dele(neop*neot*neoz),meshz(0:mpsi),&
       marker(neop*neot*neoz),adum(neop*neot*neoz),uflow(neop*neot*neoz),tp_inv,rg,dp1,delr,delp(mpsi)

  ap_inv=1.0/apart
  tp_inv=1.0/tppp(1,1)
  delr=1.0/deltar
  delp=1.0/deltap
  vp0=sqrt(2.0*tppp(1,1)*ap_inv)
  if(present(ncyclep))then
    nupp=1.88*real(pcoll)*tstep/(taupp*2.0*real(ncyclep))
  else
    nupp=1.88*real(pcoll)*tstep/taupp
  endif
  
  meshz=1.0
  if(present(z_eff))meshz=meshze
          
  vmin=1.0e-10*vp0*vp0
  dneop=real(neop)/(psi1-psi0)
  dneot=real(neot)/pi2
  dneoz=real(neoz)/pi2
 
! GC cell location, B-field, & velocity
!$omp parallel do private(m,psitmp,thetatmp,zetatmp,ip,jt,kz,isp,dpx,dp2,jst,dtx,dt2,&
!$omp& b,upara,vperp2)
  do m=1,mp
     psitmp=zpart(1,m)
     thetatmp=zpart(2,m)
     zetatmp=zpart(3,m)
     
! neoclassical cell in psi,theta, zeta
     ip=max(1,min(neop,ceiling((psitmp-psi0)*dneop)))
     jt=max(1,min(neot,ceiling(thetatmp*dneot)))           
     kz=max(1,min(neoz,ceiling(zetatmp*dneoz)))         
! GC cell
     mcell(m)=(kz-1)*neop*neot+(jt-1)*neop+ip

! radial spline grid
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
!     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline grid
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx
     
! B-field 2D spline in (psi, theta)
     b= bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
       (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
       (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

     upara=zpart(4,m)*b*qpart*ap_inv
     vperp2=zpart(6,m)*zpart(6,m)*b*2.0*ap_inv
     
! temporal storage of GC B-field, parallel & perpendicular velocity
     bgc(m)=b
     zpart(4,m)=upara
     zpart(6,m)=vperp2
  enddo

! center-of-mass coordinats; required OpenMP vector parallelization
  uflow=0.0
  marker=0.0
  do m=1,mp
     upara=zpart(4,m)
     ip=mcell(m)
     uflow(ip)=uflow(ip)+upara
     marker(ip)=marker(ip)+1.0
  enddo

! global sum
  icount=neop*neot*neoz
  call MPI_ALLREDUCE(uflow,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  uflow=adum
  call MPI_ALLREDUCE(marker,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  marker=max(1.0,adum)
  uflow=uflow/marker

! random number for ion-ion collision
  call random_number(random1(1:mp))
  call random_number(random2(1:mp))

! velocity changes due to ion-ion collision
!$omp parallel do private(m,psitmp,isp,dpx,dp2,rg,ii,dp1,pden,ptemp,ip,upara,vperp2,v,vp,zv,xv,&
!$omp& freq,k,phix,dphi,f,g,h,sp,sn,dp,dn,dpn,delu,delv2,zeff)
  do m=1,mp

! local ion density & temporature
     psitmp=zpart(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
     dp2=dpx*dpx  
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2

     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     pden=dp1*meshn(ii-1)+(1.0-dp1)*meshn(ii)
     ptemp=(dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii))*tp_inv
     zeff=dp1*meshz(ii-1)+(1.0-dp1)*meshz(ii)

! transform to center-of-mass for like-species collision
     ip=mcell(m)
     zpart(4,m)=zpart(4,m)-uflow(ip)
     upara=zpart(4,m)
     vperp2=zpart(6,m)
     v=sqrt(upara*upara+vperp2)

! velocity normalized by local temporature
     vp=vp0*sqrt(ptemp)
     zv=max(0.1,min(10.0,v/vp))
     xv=zv*zv
     
! ion-ion collision frequency
     freq=zeff*nupp*pden/(zv*xv*ptemp**1.5)

! Maxwellian integral by table look-up for range of (0.0001,10.0)
     k = min(100000,int(xv*10000.0) + 1)
     phix = (real(k)-10000.0*xv)*maxwell(k) + (10000.0*xv-real(k-1))*maxwell(k+1)
     if(xv < 0.025)phix=4.0/3.0*sqrt(xv/pi)*xv*(1.0-0.6*xv+3.0/14.0*xv*xv)
     if(xv > 10.0)phix=1.0-2.0/exp(xv)*sqrt(xv/pi)*(1.0+1.0/(2.0*xv)-1.0/(4.0*xv*xv))
     dphi = 2.0*sqrt(xv/pi)/exp(xv)

! coefficients for like-species collisions
     f=freq*2.0*phix
     g=freq*(phix-0.5*phix/xv+dphi)
     h=freq*phix/xv
     
     sp=upara*f
     sn=vperp2*(2.0*f-h-g)-2.0*upara*upara*g
     dp=upara*upara*h+vperp2*g
     dn=4.0*vperp2*v*v*v*v*g*h/dp
     dpn=2.0*vperp2*upara*(h-g)
     
! parallel and perpendicular drag and diffusion
     delu= (random1(m)-0.5)*sqrt(12.0*dp)-sp
     delv2=(random1(m)-0.5)*dpn*sqrt(12.0/dp)+(random2(m)-0.5)*sqrt(12.0*dn)-sn

! temporal storage of velocity changes
     random1(m)=delu
     random2(m)=delv2
  enddo

! momentum and energy changes due to collisions, require OpenMP vector parallelization
  dele=0.0
  delm=0.0
  do m=1,mp
     ip=mcell(m)
     delm(ip)=delm(ip)+zpart(5,m)*random1(m)
     dele(ip)=dele(ip)+zpart(5,m)*(random2(m)+random1(m)*(2.0*zpart(4,m)+random1(m)))
  enddo

! global sum
  call MPI_ALLREDUCE(delm,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  delm=adum
  call MPI_ALLREDUCE(dele,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  dele=adum
  delm=sqrt(4.5*pi)*2.0*delm/(vp0*vp0*marker)
  dele=sqrt(4.5*pi)*dele/(1.5*vp0*vp0*marker)

! local conservation of momentum and energy
!$omp parallel do private(m,upara,vperp2,psitmp,isp,dpx,dp2,rg,ii,dp1,ptemp,vp,v,zv,xv,k,phix,dphi,ip)
  do m=1,mp
! new parallel and perpendicular velocity
     upara=zpart(4,m)+random1(m)
     vperp2=max(vmin,zpart(6,m)+random2(m))

! local ion temporature
     psitmp=zpart(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
     dp2=dpx*dpx     

     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     ptemp=(dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii))*tp_inv

     vp=vp0*sqrt(ptemp)
     v=sqrt(upara*upara+vperp2)
     zv=max(0.1,min(10.0,v/vp))
     xv=zv*zv
    
     k = min(100000,int(xv*10000.0) + 1)
     phix = (k-10000.0*xv)*maxwell(k) + (10000.0*xv-k+1)*maxwell(k+1)
     if(xv .lt. 0.025)phix=4.0/3.0*sqrt(xv/pi)*xv*(1.0-0.6*xv+3.0/14.0*xv*xv)
     if(xv .gt. 10.0)phix=1.0-2.0/exp(xv)*sqrt(xv/pi)*(1.0+1.0/(2.0*xv)-1.0/(4.0*xv*xv))
     dphi = 2.0*sqrt(xv/pi)/exp(xv)
     
     ip=mcell(m)     
     zpart(5,m)=zpart(5,m)-phix/(xv*zv)*upara*delm(ip)-(phix-dphi)/zv*dele(ip)
     zpart(4,m)=(upara+uflow(ip))*apart/(qpart*bgc(m))
     zpart(6,m)=sqrt(0.5*vperp2*apart/bgc(m))
  enddo

end subroutine collision_fokkerplanck
