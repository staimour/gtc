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

subroutine push(species_name,icycle,irksub,ihybrid)
  use global_parameters,only: nparami,nparame,nparamf,nparamfe,ncyclee,ncyclefe
  use particle_array
  use global_parameters,only: mype
  implicit none

  integer ierror
  character(*),intent(in) :: species_name
  integer,intent(in),optional :: icycle,irksub,ihybrid

  interface
    subroutine gkPushParticle(zpart,zpart0,wppart,wtpart0,wtpart1,&
        data1d,jtpart0,jtpart1,diagpart,wzpart,rdtem,kapan,kapat,meshn,mesht,&
        pflux,qpart,apart,pload,ngyro,mpmax,mp,nparam,ncyclep,icycle,irksub,&
        mp1,zpart1,ihybrid,phit,dnt)
      use precision
      use particle_array,only: mpdiag
      implicit none

      !declaration of the dummy arguments
      integer pload,ngyro,mpmax,mp,nparam
      integer,optional :: ncyclep,icycle,irksub,ihybrid,mp1
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart
      real(wp),dimension(0:) :: kapan,kapat,meshn,mesht,pflux,rdtem
      real(wp),dimension(mpdiag) :: diagpart
      real(wp),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: data1d
      real(wp),dimension(:,:),optional :: zpart1
      real(wp),dimension(0:,:),optional :: phit,dnt
    end subroutine gkPushParticle

    subroutine hybridPushParticle(zpart,zpart0,wppart,wtpart0,wtpart1,&
        data1d,jtpart0,jtpart1,diagpart,wzpart,rdtem,kapan,kapat,meshn,mesht,&
        pflux,qpart,apart,pload,ngyro,mpmax,mp,nparam,ncyclep,icycle,irksub,&
        mp1,zpart1,ihybrid,phit,dnt)
      use precision
      use particle_array,only: mpdiag
      implicit none

      integer pload,ngyro,mpmax,mp,nparam
      integer,optional :: ncyclep,icycle,irksub,ihybrid,mp1
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart
      real(wp),dimension(0:) :: kapan,kapat,meshn,mesht,pflux,rdtem
      real(wp),dimension(mpdiag) :: diagpart
      real(wp),dimension(:,:) :: zpart0,wppart,wtpart0,wtpart1
      real(wp),dimension(:,:) :: zpart
      real(wp),dimension(0:,:) :: data1d
      real(wp),dimension(:,:),optional :: zpart1
      real(wp),dimension(0:,:),optional :: phit,dnt
    end subroutine hybridPushParticle
  end interface

  if(species_name=="thermal-ion")then
    call gkPushParticle(zion,zion0,wpion,wtion0,wtion1,data1di,&
      jtion0,jtion1,diagion,wzion,rdtemi,kapani,kapati,meshni,meshti,pfluxi,&
      qion,aion,iload,ngyroi,mimax,mi,nparami)
  elseif(species_name=="thermal-electron")then
    call hybridPushParticle(zelectron,zelectron0,wpelectron,wtelectron0,&
      wtelectron1,data1de,jtelectron0,jtelectron1,diagelectron,wzelectron,&
      rdteme,kapane,kapate,meshne,meshte,pfluxe,qelectron,aelectron,eload,&
      ngyroe,memax,me,nparame,ncyclee,icycle,irksub,me1,zelectron1,ihybrid,&
      phit,dnet)
  elseif(species_name=="fast-ion")then
    call gkPushParticle(zfast,zfast0,wpfast,wtfast0,wtfast1,data1df,&
      jtfast0,jtfast1,diagfast,wzfast,rdtemf,kapanf,kapatf,meshnf,meshtf,&
      pfluxf,qfast,afast,fload,ngyrof,mfmax,mf,nparamf)
#ifdef _FRC
  elseif(species_name=="fast-electron")then
    call gkPushParticle(zfaste,zfaste0,wpfaste,wtfaste0,wtfaste1,&
     data1dfe,jtfaste0,jtfaste1,diagfaste,wzfaste,rdtemfe,kapanfe,kapatfe,&
     meshnfe,meshtfe,pfluxfe,qfaste,afaste,feload,ngyrofe,mfemax,mfe,nparamfe)
#else
  elseif(species_name=="fast-electron")then
    call gkPushParticle(zfaste,zfaste0,wpfaste,wtfaste0,wtfaste1,&
     data1dfe,jtfaste0,jtfaste1,diagfaste,wzfaste,rdtemfe,kapanfe,kapatfe,&
     meshnfe,meshtfe,pfluxfe,qfaste,afaste,feload,ngyrofe,mfemax,mfe,&
     nparamfe,ncyclefe,icycle,irksub,mfe1,zfaste1)
#endif
  else
    write(*,*)'push.F90: wrong choice: ',species_name
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine push
