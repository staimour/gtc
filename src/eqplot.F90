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

subroutine eqplot
  use global_parameters
  use field_array
  use equilibrium
  implicit none
  
  integer,parameter :: mskip=4
  integer i,j,ieq,nplot,nrad
  real delpsi,dum,pdum,tdum,deltheta,datap(lsp),data1d(lsp),&
    datax(mpsi/mskip+1,lst),dataz(mpsi/mskip+1,lst),data2d(mpsi/mskip+1,lst,5)
  real,external :: splcos,splsin,spb,spx,spz,spjac,currenti,zeta2phi,delb,&
    spzeff,sptrot,sper,spq,spgpsi,sprpsi,sprgpsi,sptorpsi,&
    sppsitor,sppsirg,sppressure,dq_dp,spden,den_dp,sptem,tem_dp
       

  ieq=120 !open equilibiurm plot output file
  open(ieq,file='equilibrium.out',status='replace')
101 format(i6)
102 format(e16.8)

!# of 1D radial plots
  nplot=29
!  nrad=mpsi+1
  nrad=lsp
  write(ieq,101)nplot,nrad

!0: radial axis using poloidal flux function
!  delpsi=(psimesh(mpsi)-psimesh(0))/real(nrad-1)
  delpsi=psiw/real(nrad-1)
  do i=1,nrad
!     datap(i)=psimesh(0)+delpsi*real(i-1)
     datap(i)=delpsi*real(i-1)
  enddo
  write(ieq,102)datap 
  
!1: sqaure-root of normalized toroidal flux function
  data1d=0.0
  write(ieq,102)data1d 

!2: minor radius
  data1d=0.0
  write(ieq,102)data1d 

!3: major radius
  data1d=0.0
  write(ieq,102)data1d 

#ifdef _FRC
!4: Te
  do i=1,nrad
     data1d(i)=sptem(datap(i),tfepp)/(rho0*rho0)
  enddo
  write(ieq,102)data1d 

!5: -d(ln(Te))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-tem_dp(dum,tfepp)/sptem(dum,tfepp)
  enddo
  write(ieq,102)data1d 

!6: ne
  do i=1,nrad
     data1d(i)=spden(datap(i),nfepp)
  enddo
  write(ieq,102)data1d 

!7: -d(ln(ne))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-den_dp(dum,nfepp)/spden(dum,nfepp)
  enddo
  write(ieq,102)data1d 
#else
!4: Te
  do i=1,nrad
     data1d(i)=sptem(datap(i),tepp)/tepp(1,1)
  enddo
  write(ieq,102)data1d 

!5: -d(ln(Te))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-tem_dp(dum,tepp)/sptem(dum,tepp)
  enddo
  write(ieq,102)data1d 

!6: ne
  do i=1,nrad
     data1d(i)=spden(datap(i),nepp)/nepp(1,1)
  enddo
  write(ieq,102)data1d 

!7: -d(ln(ne))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-den_dp(dum,nepp)/spden(dum,nepp)
  enddo
  write(ieq,102)data1d 
#endif
!8: ti
  do i=1,nrad
     data1d(i)=sptem(datap(i),tipp)/tepp(1,1)
  enddo
  write(ieq,102)data1d 

!9: -d(ln(ti))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-tem_dp(dum,tipp)/sptem(dum,tipp)
  enddo
  write(ieq,102)data1d

!10: ni
  do i=1,nrad
     data1d(i)=spden(datap(i),nipp)/nepp(1,1)
  enddo
  write(ieq,102)data1d 

!11: -d(ln(ni))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-den_dp(dum,nipp)/spden(dum,nipp)
  enddo
  write(ieq,102)data1d 

#ifndef _FRC
!12: tf
  do i=1,nrad
     data1d(i)=sptem(datap(i),tfpp)/tepp(1,1)
  enddo
  write(ieq,102)data1d 

!13: -d(ln(tf))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-tem_dp(dum,tfpp)/sptem(dum,tfpp)
  enddo
  write(ieq,102)data1d

!14: nf
  do i=1,nrad
     data1d(i)=spden(datap(i),nfpp)/nepp(1,1)
  enddo
  write(ieq,102)data1d 

!15: -d(ln(nf))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-den_dp(dum,nfpp)/spden(dum,nfpp)
  enddo
  write(ieq,102)data1d 
#else
!12: tf
  do i=1,nrad
     data1d(i)=sptem(datap(i),tepp)/(rho0*rho0)
  enddo
  write(ieq,102)data1d 

!13: -d(ln(tf))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-tem_dp(dum,tepp)/sptem(dum,tepp)
  enddo
  write(ieq,102)data1d

!14: nf
  do i=1,nrad
     data1d(i)=spden(datap(i),nepp)
  enddo
  write(ieq,102)data1d 

!15: -d(ln(nf))/dr
  do i=1,nrad
     dum=datap(i)
     data1d(i)=-den_dp(dum,nepp)/spden(dum,nepp)
  enddo
  write(ieq,102)data1d 
#endif
!16: zeff
  do i=1,nrad
     data1d(i)=spzeff(datap(i))
  enddo
  write(ieq,102)data1d 

!17: toroidal rotation
  do i=1,nrad
     data1d(i)=sptrot(datap(i))
  enddo
  write(ieq,102)data1d 

!18: radial electric field
  do i=1,nrad
     data1d(i)=sper(datap(i))
  enddo
  write(ieq,102)data1d 

!19: q profile
  do i=1,nrad
     data1d(i)=spq(datap(i))
  enddo
  write(ieq,102)data1d 

!20: d(ln(q))/dpsi
  do i=1,nrad
     dum=datap(i)
     data1d(i)=dq_dp(dum)/spq(dum)
  enddo
  write(ieq,102)data1d

!21: gcurrent profile
  do i=1,nrad
     data1d(i)=spgpsi(datap(i))
  enddo
  write(ieq,102)data1d 

!22: pressure profile
  do i=1,nrad
     data1d(i)=sppressure(datap(i))
  enddo
  write(ieq,102)data1d 

!23: minor radius
  do i=1,nrad
     data1d(i)=sprpsi(datap(i))
  enddo
  write(ieq,102)data1d 

!24: toroidal flux
  do i=1,nrad
     data1d(i)=sptorpsi(datap(i))
  enddo
  write(ieq,102)data1d

!25: radial grid: rgpsi
  do i=1,nrad
     data1d(i)=sprgpsi(datap(i))
  enddo
  write(ieq,102)data1d

!26: inverse of spline torpsi: psitor
  tdum=0.0
  pdum=sptorpsi(psiw)
  dum=(pdum-tdum)/real(nrad-1)
  do i=1,nrad
     pdum=tdum+dum*real(i-1)
     data1d(i)=sppsitor(pdum)
  enddo
  write(ieq,102)data1d

!27: inverse of spline rgpsi: psirg
  tdum=0.0
  pdum=sprgpsi(psiw)
  dum=(pdum-tdum)/real(nrad-1)
  do i=1,nrad
     pdum=tdum+dum*real(i-1)
     data1d(i)=sppsirg(pdum)
  enddo
  write(ieq,102)data1d

!28: error of spline cos in [0, pi/2]
  do i=1,nrad
     dum=0.25*spdtheta*real(i-1)
     data1d(i)=splcos(dum)-cos(dum)
  enddo
  write(ieq,102)data1d

!29: error of spline sin in [0, pi/2]
  do i=1,nrad
     dum=0.25*spdtheta*real(i-1)
     data1d(i)=splsin(dum)-sin(dum)
  enddo
  write(ieq,102)data1d

!# of 2D contour plots on poloidal plane
  nplot=5
  write(ieq,101)nplot,mpsi/mskip+1,lst
  deltheta=2.0*pi/real(lst-1)

  do i=1,mpsi/mskip+1
     pdum=psimesh((i-1)*mskip)
     do j=1,lst
        tdum=deltheta*real(j-1)
        datax(i,j)=spx(pdum,tdum)
        dataz(i,j)=spz(pdum,tdum)
        data2d(i,j,1)=spb(pdum,tdum)
        data2d(i,j,2)=spjac(pdum,tdum)
        data2d(i,j,3)=currenti(pdum,tdum)
        data2d(i,j,4)=zeta2phi(pdum,tdum)
        data2d(i,j,5)=delb(pdum,tdum)
     enddo
  enddo

!0-1: mesh points on (X,Z)
  write(ieq,102)datax,dataz 

!2: b-field
  write(ieq,102)data2d(:,:,1)

!3: Jacobian
  write(ieq,102)data2d(:,:,2)

!4: icurrent
  write(ieq,102)data2d(:,:,3)

!5: zeta2phi
  write(ieq,102)data2d(:,:,4)

!6: delb
  write(ieq,102)data2d(:,:,5)

close(ieq)
end subroutine eqplot
