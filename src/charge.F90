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

subroutine charge(species_name)
  use global_parameters,only: ihybrid,nhybrid
  use particle_array
  use field_array,only: sfluidne
  implicit none

  interface
    subroutine gkChargeParticle(zpart,wppart,wtpart0,wtpart1,density,&
        flow,jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,&
        qpart,apart,pload,ngyro,mp,ihybrid,nhybrid,&
        pressurepara,pressureperp,sfluidn,dnsave,&
        switch)
      use precision
      implicit none

      integer pload,ngyro,mp
      integer,optional :: ihybrid,nhybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart,marker,markert
      real(wp),dimension(0:) :: zonal,zonalc,pmark
      real(wp),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: density,flow
      real(wp),dimension(:,:) :: zpart
      real(wp),dimension(:,:),optional :: pressurepara,pressureperp,sfluidn
      real(wp),dimension(:,:,:),optional :: dnsave
      character(*),intent(in),optional :: switch
    end subroutine gkChargeParticle

    subroutine hybridChargeParticle(zpart,wppart,wtpart0,wtpart1,density,&
        flow,jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,&
        qpart,apart,pload,ngyro,mp,ihybrid,nhybrid,&
        pressurepara,pressureperp,sfluidn,dnsave,&
        switch)
      use precision
      implicit none

      integer pload,ngyro,mp
      integer,optional :: ihybrid,nhybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart,marker,markert
      real(wp),dimension(0:) :: zonal,zonalc,pmark
      real(wp),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: density,flow
      real(wp),dimension(:,:) :: zpart
      real(wp),dimension(:,:),optional :: pressurepara,pressureperp,sfluidn
      real(wp),dimension(:,:,:),optional :: dnsave
      character(*),intent(in),optional :: switch
    end subroutine hybridChargeParticle
  end interface

  integer ierror
  character(*),intent(in) :: species_name

  if(species_name=='thermal-ion')then
    !ideal MHD
    if(iload==0)then
      densityi=0.0
      flowi=0.0
    else
      call gkChargeParticle(zion,wpion,wtion0,wtion1,densityi,flowi,jtion0,&
        jtion1,wzion,zonali,zonalci,markeri,markerit,pmarki,qion,aion,&
        iload,ngyroi,mi,switch="w/o density modification")
    endif
  elseif(species_name=='fast-ion')then
    call gkChargeParticle(zfast,wpfast,wtfast0,wtfast1,densityf,flowf,&
      jtfast0,jtfast1,wzfast,zonalf,zonalcf,markerf,markerft,pmarkf,&
      qfast,afast,fload,ngyrof,mf)
  elseif(species_name=='fast-electron')then
    call gkChargeParticle(zfaste,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,&
      jtfaste0,jtfaste1,wzfaste,zonalfe,zonalcfe,markerfe,markerfet,pmarkfe,&
      qfaste,afaste,feload,ngyrofe,mfe)
  elseif(species_name=='thermal-electron')then
    call hybridChargeParticle(zelectron,wpelectron,wtelectron0,wtelectron1,&
      densitye,flowe,jtelectron0,jtelectron1,wzelectron,zonale,zonalce,&
      markere,markeret,pmarke,qelectron,aelectron,eload,ngyroe,me,&
      ihybrid=ihybrid,nhybrid=nhybrid,pressurepara=pressureepara,&
      pressureperp=pressureeperp,sfluidn=sfluidne,dnsave=dnesave)
  else
    write(*,*)'push.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine charge
