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
! ==============================================================================

! function of plasma profile
! electron density
real function spden(pdum,denpp)
  use equilibrium, only: lsp, spdpsi  
  implicit none
  real pdum
  real,dimension(3,lsp) :: denpp
  real,external :: spline0

  spden=spline0(pdum,lsp,spdpsi,denpp)
end function spden

! d_densitye/d_psi
real function den_dp(pdum,denpp)
  use equilibrium, only: lsp, spdpsi  
  implicit none
  real pdum
  real,dimension(3,lsp) :: denpp
  real,external :: dspline0

  den_dp=dspline0(pdum,lsp,spdpsi,denpp)
end function den_dp

! electron temperature
real function sptem(pdum,tpp)
  use equilibrium, only: lsp, spdpsi  
  implicit none
  real pdum
  real,dimension(3,lsp) :: tpp
  real,external :: spline0

  sptem=spline0(pdum,lsp,spdpsi,tpp)
end function sptem

! d_etemperature/d_psi
real function tem_dp(pdum,tpp)
  use equilibrium, only: lsp, spdpsi  
  implicit none
  real pdum
  real,dimension(3,lsp) :: tpp
  real,external :: dspline0

  tem_dp=dspline0(pdum,lsp,spdpsi,tpp)
end function tem_dp

! Z_eff
real function spzeff(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  spzeff=spline0(pdum,lsp,spdpsi,zepp)
end function spzeff

! toroidal rotation
real function sptrot(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  sptrot=spline0(pdum,lsp,spdpsi,ropp)
end function sptrot

! E_r
real function sper(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  sper=spline0(pdum,lsp,spdpsi,erpp)
end function sper

! spline cos function
real function splcos(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  splcos=spline0(pdum,lst,spdtheta,spcos)
end function splcos

! spline sin function
real function splsin(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  splsin=spline0(pdum,lst,spdtheta,spsin)
end function splsin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function from equilibrium data
! safety factor
real function spq(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  spq=spline0(pdum,lsp,spdpsi,qpsi)
end function spq

! dq_dpsi
real function dq_dp(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: dspline0

  dq_dp=dspline0(pdum,lsp,spdpsi,qpsi)
end function dq_dp

! g and dg/dp
real function spgpsi(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  spgpsi=spline0(pdum,lsp,spdpsi,gpsi)
end function spgpsi

real function spgpsi_dp(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: dspline0

  spgpsi_dp=dspline0(pdum,lsp,spdpsi,gpsi)
end function spgpsi_dp

! I and dI/dp
real function spcpsi(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  spcpsi=spline0(pdum,lsp,spdpsi,cpsi)
end function spcpsi

real function spicpsi_dp(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: dspline0

  spicpsi_dp=dspline0(pdum,lsp,spdpsi,cpsi)
end function spicpsi_dp

! pressure
real function sppressure(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  sppressure=spline0(pdum,lsp,spdpsi,ppsi)
end function sppressure

! poloidal to toroidal flux
real function sptorpsi(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  sptorpsi=spline0(pdum,lsp,spdpsi,torpsi)
end function sptorpsi

! d(toroidal flux)/d(poloidal flux)
real function dtorpsi(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: dspline0

  dtorpsi=dspline0(pdum,lsp,spdpsi,torpsi)
end function dtorpsi

! toroidal to poloidal flux
real function sppsitor(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline0

  sppsitor=spline0(pdum,lsp,spdtor,psitor)
end function sppsitor

! minor radius
real function sprpsi(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline1

  sprpsi=spline1(pdum,lsp,spdpsi,rpsi)
end function sprpsi

! dr/dpsi
real function dr_dp(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: dspline1

  dr_dp=dspline1(pdum,lsp,spdpsi,rpsi)
end function dr_dp

! psi to radial grid
real function sprgpsi(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline1

  sprgpsi=spline1(pdum,lsp,spdpsi,rgpsi)
end function sprgpsi

! radial grid to psi
real function sppsirg(pdum)
  use equilibrium
  implicit none
  real pdum
  real,external :: spline1

  sppsirg=spline1(pdum,lsp,spdrg,psirg)
end function sppsirg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! magnetic field amplitude
real function spb(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  spb=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,bsp)
!  spb=1.0 !cylindrical limit

end function spb

! db_field/dpsi
real function dbdp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dbdp=spline2d(1,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,bsp)
end function dbdp

! db_field/dtheta
real function dbdt(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dbdt=spline2d(2,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,bsp)
end function dbdt

! transform from (psi,theta) to (X,Z)
real function spx(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  spx=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,xsp)
end function spx

! transform from (psi,theta) to (X,Z)
real function spz(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  spz=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,zsp)
end function spz

!transformation dR/dpsi
real function dxdp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dxdp=spline2d(1,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,xsp)
end function dxdp

!dz/dpsi
real function dzdp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dzdp=spline2d(1,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,zsp)
end function dzdp

!dx/dtheta
real function dxdt(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dxdt=spline2d(2,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,xsp)
end function dxdt

!dz/dtheta
real function dzdt(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dzdt=spline2d(2,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,zsp)
end function dzdt
!!!!!!!!!!!!!!!!!!!--3D-case--!!!!!!!!!!!!!!!!!!!!!!!!!!!
!transformation d\Phi/dpsi
real function dfdp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dfdp=spline2d(1,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,fsp)
end function dfdp

!transformation d\Phi/dtheta
real function dfdt(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dfdt=spline2d(2,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,fsp)
end function dfdt

!transformation d\Phi/dzeta
real function dfdz(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dfdz=1.0+spline2d(3,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,fsp) !!! 1.0 since \Phi is determined as delta\Phi
end function dfdz

!dx/dzeta
real function dxdz(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dxdz=spline2d(3,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,xsp)
end function dxdz

!dz/dzeta
real function dzdz(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dzdz=spline2d(3,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,zsp)
end function dzdz

! Jacobian
real function spjac(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  spjac=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,gsp)
end function spjac

! toroidal current I
real function currenti(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  currenti=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,rd)
end function currenti

real function didp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  didp=spline2d(1,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,rd)
end function didp

! difference between magnetic angle zeta and cylindrical angle phi
real function zeta2phi(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  zeta2phi=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,nu)
end function zeta2phi

! delta in B-field contravariant representation
real function delb(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  delb=spline2d(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,dl)
end function delb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1D spline for radial profiles
real function spline0(pdum,nsp,delx,y)
  implicit none
  integer nsp,i
  real pdum,y(3,nsp),delx,dpx

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  spline0=y(1,i)+dpx*y(2,i)+dpx*dpx*y(3,i)

end function spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of 1D spline function
real function dspline0(pdum,nsp,delx,y)
  implicit none
  integer nsp,i
  real pdum,y(3,nsp),delx,dpx

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  dspline0=y(2,i)+2.0*dpx*y(3,i)

end function dspline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1D spline with first point being linear function y=sqrt(x)
real function spline1(pdum,nsp,delx,y)
  implicit none
  integer nsp,i
  real pdum,y(3,nsp),delx,dpx,dp2

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
! expand y(x) using sprt(x) near x=0
  if(i==1)dpx=sqrt(dpx)
  dp2=dpx*dpx
  spline1=y(1,i)+dpx*y(2,i)+dp2*y(3,i)

end function spline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of 1D spline with first point being linear function y=sqrt(x)
real function dspline1(pdum,nsp,delx,y)
  implicit none
  integer nsp,i
  real pdum,y(3,nsp),delx,dpx

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  if(i==1)dpx=sqrt(dpx)

  if(i==1)then
     dspline1=0.5*y(2,i)/dpx+y(3,i) ! y(x)=y1+y2*sprt(x)+y3*x near x=0
  else
     dspline1=y(2,i)+2.0*dpx*y(3,i) ! y(x)=y1+y2*x+y3*x*x otherwise
  endif

end function dspline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2D spline on poloidal plane
real function spline2d(iflag,sd,pdum,tdum,nsp,nst,delp,delt,f)
  implicit none
  integer iflag,nsp,nst,i,j,is,sd
  real pdum,tdum,f(sd,nsp,nst),dpx,dp2,dtx,dt2,delp,delt,dx(sd)

  i=max(1,min(nsp-1,ceiling(pdum/delp)))
  dpx=pdum-delp*real(i-1)
! remove sigularity near psi=1 due to functional form of sqrt(psi)
! expand f(x,y) using sprt(x) near x=0; i.e., subroutine spline1 in 1D case
  if(i==1)dpx=sqrt(dpx)
  dp2=dpx*dpx

  j=max(1,min(nst-1,ceiling(tdum/delt)))
  dtx=tdum-delt*real(j-1)
  dt2=dtx*dtx

  dx=0.0
  dx(1)=1.0
  dx(2)=dpx
  dx(3)=dp2
  dx(4:6)=dx(1:3)*dtx
  dx(7:9)=dx(1:3)*dt2

  spline2d=0.0
  if(iflag==0)then !2D spline value
     spline2d=f(1,i,j)    +f(2,i,j)*dpx    +f(3,i,j)*dp2 &
          +f(4,i,j)*dtx+f(5,i,j)*dtx*dpx+f(6,i,j)*dtx*dp2 &
          +f(7,i,j)*dt2+f(8,i,j)*dt2*dpx+f(9,i,j)*dt2*dp2
  elseif(iflag==1)then !derivative with respect to psi
     if(i==1)then
        spline2d=0.5*(f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2)/dpx+f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2
     else
        spline2d=f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2+2.0*dpx*(f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2)
     endif
  elseif(iflag==2) then !derivative with respect to theta
     spline2d=f(4,i,j)+f(5,i,j)*dpx+f(6,i,j)*dp2+2.0*dtx*(f(7,i,j)+f(8,i,j)*dpx+f(9,i,j)*dp2)
  elseif(iflag==3) then !derivative with respect to zeta
     if(sd==27)then
       dx(10:18)=dx(1:9)
       do is = 10, 18
          spline2d=spline2d+f(is,i,j)*dx(is)
       enddo
     else
       spline2d=0.0
     endif
  endif
end function spline2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! construct 1D spline with on x=[0,xmax], grid x_i = (i-1)*delx, delx=xmax/(nsp-1)
! in domain i, y = y(1,i) + y(2,i)*delx + y(3,i)*delx**2
subroutine construct_spline0(iflag,nsp,delx,y)
  implicit none
  integer i,nsp,ipp,iflag
  real delx,y(3,nsp)

! first point
  if(iflag==0)then
! iflag=0: first point being y=y1+y2*x+y3*x*x
     y(2,1)=(4.0*y(1,2)-y(1,3)-3.0*y(1,1))/(2.0*delx)
     y(3,1)=(y(1,2)-y(1,1)-y(2,1)*delx)/(delx*delx)

  elseif(iflag==1)then
! iflag=1: first point being linear function y=y_1+y2*x
     y(2,1)=(y(1,2)-y(1,1))/delx
     y(3,1)=0.0

  elseif(iflag==2)then
! iflag=2: first point being quadratic function y=y1+y3*x*x
     y(2,1)=0.0
     y(3,1)=(y(1,2)-y(1,1))/(delx*delx)
  endif

  do i=2,nsp-2
     ipp=min(i+2,nsp)
     y(2,i)=-y(2,i-1)+2.0*(y(1,i)-y(1,i-1))/delx

! smooth f1
     y(1,i+1)=0.5*delx*y(2,i)+0.25*y(1,ipp)+0.75*y(1,i)
  enddo

  y(2,nsp-1)=-y(2,nsp-2)+2.0*(y(1,nsp-1)-y(1,nsp-2))/delx
  y(2,nsp)=-y(2,nsp-1)+2.0*(y(1,nsp)-y(1,nsp-1))/delx

  do i=2,nsp-1
     y(3,i)=(y(2,i+1)-y(2,i))/(2.0*delx)
  enddo

! last point is not used;
  y(3,nsp)=0.0

end subroutine construct_spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inversion of spline0 y(x) to x(y)
subroutine invert_spline0(iflag,nsp,delx,dely,y,x)
  implicit none
  integer i,nsp,j,iflag
  real delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

! first point given by input
  x(1,1)=0.0
! other points
  do i=2,nsp-1
! y grid
     ydum=dely*real(i-1)
! search x grid for ydum
     j=1
     do while (ydum>y(1,j+1))
        j=j+1
     enddo

! x(1,i)=j grid location + distance from j grid
     y0=y(1,j)
     y1=y(2,j)
     y2=y(3,j)
     if (abs(y2)>0.000001) then
       x(1,i)=delx*real(j-1)+(sqrt(y1*y1+4.0*y2*(ydum-y0))-y1)/(2.0*y2)
     else
       x(1,i)=delx*real(j-1)+(ydum-y0)/y1
     endif

  enddo

! last point
  x(1,nsp)=x(1,1)+delx*real(nsp-1)

! spline fit x function
  if(iflag==0)call construct_spline0(0,nsp,dely,x)

end subroutine invert_spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! spline1 for cases with first point being function y=y1+y2*sqrt(x)+y3*x
subroutine construct_spline1(nsp,delx,f)
  implicit none
  integer i,nsp,ipp
  real delx,f(3,nsp)

! first point
  f(2,1)=(2.0*f(1,2)-f(1,3)-f(1,1))/((2.0-sqrt(2.0))*sqrt(delx))
  f(3,1)=(f(1,2)-f(1,1)-f(2,1)*sqrt(delx))/delx

! second point
  f(2,2)=0.5*f(2,1)/sqrt(delx)+f(3,1)
  f(3,2)=(f(1,3)-f(1,2)-delx*f(2,2))/(delx*delx)

  do i=3,nsp-2
     ipp=min(i+2,nsp)
     f(2,i)=-f(2,i-1)+2.0*(f(1,i)-f(1,i-1))/delx

! smooth f1
     f(1,i+1)=0.5*delx*f(2,i)+0.25*f(1,ipp)+0.75*f(1,i)
  enddo

  f(2,nsp-1)=-f(2,nsp-2)+2.0*(f(1,nsp-1)-f(1,nsp-2))/delx
  f(2,nsp)=-f(2,nsp-1)+2.0*(f(1,nsp)-f(1,nsp-1))/delx

  do i=3,nsp-1  !!change start from 1 to 3
     f(3,i)=(f(2,i+1)-f(2,i))/(2.0*delx)
  enddo

! last point is not used;
  f(3,nsp)=0.0

end subroutine construct_spline1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inversion of spline1 y(x) to x(y)
subroutine invert_spline1(nsp,delx,dely,y,x)
  implicit none
  integer i,nsp,j
  real delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

! first point is given by inputs
  x(1,1)=0.0
! other points
  do i=2,nsp-1
! y grid
     ydum=dely*real(i-1)
! search x grid for ydum
     j=1
     do while (ydum>y(1,j+1))
        j=j+1
     enddo

! x(1,i)=j grid location + distance from j grid
     y0=y(1,j)
     y1=y(2,j)
     y2=y(3,j)

!     if(j==1)then
!        x(1,i)=((ydum-y0)/y(2,1))**2 ! first point uses sqrt(x)
!     else
!        x(1,i)=delx*real(j-1)+(sqrt(y1*y1+4.0*y2*(ydum-y0))-y1)/(2.0*y2)
!     endif
     if (abs(y2)>0.000001) then
       x(1,i)=delx*real(j-1)+(sqrt(y1*y1+4.0*y2*(ydum-y0))-y1)/(2.0*y2)
     else
       x(1,i)=delx*real(j-1)+(ydum-y0)/y1
     endif

     if (j==1) x(1,i)=x(1,i)*x(1,i)

  enddo
! last point
  x(1,nsp)=delx*real(nsp-1)

! spline fit x function
! call spline with first point being ~ sqrt(r)
  call construct_spline1(nsp,dely,x)

end subroutine invert_spline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine construct_spline2d(nx,ny,delx,dely,y)
  implicit none
  integer i,j,s,nx,ny,ipp
  real delx,dely,y(9,nx,ny),ddum1(3,nx),ddum2(3,ny)

  ddum1=0.0
  ddum2=0.0

  do j = 1, ny
     ddum1(1,:)=y(1,:,j)
     call construct_spline1(nx,delx,ddum1)
     y(1,:,j)=ddum1(1,:)
     y(2,:,j)=ddum1(2,:)
     y(3,:,j)=ddum1(3,:)
  enddo
  do i = 2, nx
     do s = 1, 3
        ddum2(1,:)=y(s,i,:)
        call construct_spline0(0,ny,dely,ddum2)
        y(s,i,:)=ddum2(1,:)
        y(s+3,i,:)=ddum2(2,:)
        y(s+6,i,:)=ddum2(3,:)
     enddo
  enddo
end subroutine construct_spline2d
!!!!!!!!!!!!!
subroutine spline3d(nx,ny,nz,delx,dely,ndim,ntor,ycn,ysn,btemp,sgn)
  use precision
  use global_parameters
  implicit none

  integer i,j,k,n,s

  integer,intent(in) :: nx,ny,nz,ndim,sgn
  integer,intent(in) :: ntor(ndim)
  real(wp),intent(in) :: delx,dely,ycn(nx,ny,ndim),ysn(nx,ny,ndim)
  real(wp),intent(out) :: btemp(27,nx,ny,nz)
  real(wp) :: delz,zdum,dum,ddum(9,nx,ny),ddum3(3,nz+1)

  delz=2.0*pi/real(nz)
  ddum=0.0
  ddum3=0.0
  btemp=0.0

! for each toroidal grid the value is calculated as sum over cos and sin harmonics
! negative sign in zeta due to VMEC convention, sgn flipps directions of theta and zeta
! to ensure poloidal flux is positive
  if(sgn>0)then
    do k = 1, nz
      zdum= -real(k-1)*delz 
      btemp(1,:,:,k)=0.0
      do i = 1, nx
          do j = 1, ny 
             dum=0.0
             do n = 1, ndim
                dum=dum+ycn(i,j,n)*cos(real(ntor(n))*zdum)+ysn(i,j,n)*sin(real(ntor(n))*zdum)
             enddo
             btemp(1,i,j,k)=dum
          enddo
      enddo
      ddum(1,:,:)=btemp(1,:,:,k)
      call construct_spline2d(nx,ny,delx,dely,ddum)
      btemp(1:9,:,:,k)=ddum(1:9,:,:)
    enddo
  else !VMEC original z-direction
    do k = 1, nz
      zdum= real(k-1)*delz ! changed sign of zeta
      btemp(1,:,:,k)=0.0
      do i = 1, nx
          do j = 1, ny ! change order for sgn<0
             dum=0.0
             do n = 1, ndim
                dum=dum+ycn(i,ny-j+1,n)*cos(real(ntor(n))*zdum)+ysn(i,ny-j+1,n)*sin(real(ntor(n))*zdum)
             enddo
             btemp(1,i,j,k)=dum 
          enddo
      enddo
      ddum(1,:,:)=btemp(1,:,:,k)
      call construct_spline2d(nx,ny,delx,dely,ddum)
      btemp(1:9,:,:,k)=ddum(1:9,:,:)
    enddo
  endif

  do i = 1, nx
    do j = 1, ny
        do s = 1, 9
            ddum3(1,:)=btemp(s,i,j,:)
            ddum3(1,nz+1)=btemp(s,i,j,1) ! extra point to ensure periodicity in spline construction
            call construct_spline0(0,nz+1,delz,ddum3)
            btemp(s,i,j,:)=ddum3(1,1:nz)
            btemp(s+9,i,j,:)=ddum3(2,1:nz)
            btemp(s+18,i,j,:)=ddum3(3,1:nz)
        enddo
    enddo
  enddo
end subroutine spline3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cross(vec1,vec2,vec3)
  implicit none
  real:: vec1(3),vec2(3),vec3(3)
  
  vec3(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
  vec3(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
  vec3(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
  
end subroutine cross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function dots(vec1,vec2,n)
  implicit none
  integer i,n
  real ::temp,vec1(n),vec2(n)
  temp=0.0
  do i=1,n
    temp=temp+vec1(i)*vec2(i)
  enddo
  dots=temp
end function dots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
