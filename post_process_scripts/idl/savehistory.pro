pro savehistory

filedirectories=["m=5_want445"];["omni_m=10_950","omni_m=10_850","omni_m=10_700","omni_m=10_800"];
numfd=n_elements(filedirectories);
for nfi=0,numfd-1 do begin
  cd,'/scratch/scratchdirs/clau/V0124/Drift/'+filedirectories[nfi]
  fname="history.out"
  
  ihistory=11
  openr, ihistory, fname
  
  ; # of time steps
  ndstep=1L
  ; # of species: ion, electron, EP, impuries
  nspecies=1L
  ; # of quantities per species: density,entropy,momentum,energy, fluxes
  mpdiag=1L
  ; # of field variables: phi, a_par, fluidne
  nfield=1L
  ; # of modes per field: (n,m)
  modes=1L
  ; # of quantities per field: rms, single spatial point, zonal components
  mfdiag=1L
  ; # time step size
  tstep=0.0
  
  readf,ihistory, ndstep,nspecies,mpdiag,nfield,modes,mfdiag,tstep
  print, ndstep,nspecies,mpdiag,nfield,modes,mfdiag,tstep
  
  ;        print,'total time steps=',ndstep
  ;        print,'read # of last time step'
  ;        new=''
  ;        read,new
  ;        if new ne '' then ndstep=0+new
  partdata=fltarr(mpdiag,nspecies)
  fieldtime=fltarr(mfdiag,nfield)
  fieldmode=fltarr(2,modes,nfield)
  
  parthist=fltarr(ndstep,mpdiag,nspecies)
  fieldhist=fltarr(ndstep,mfdiag,nfield)
  modehist=fltarr(ndstep,2,modes,nfield)
  
  
  on_ioerror,ioerr
  i=0L
  while i lt ndstep do begin
     ;for i=0,ndstep-1 do begin
        readf,ihistory,partdata,fieldtime,fieldmode
        parthist(i,*,*)=partdata
        fieldhist(i,*,*)=fieldtime
        modehist(i,*,*,*)=fieldmode
    
    ; print,i+1,modehist(i,0,4,2),modehist(i,1,4,2)
     ;endfor
    i=i+1
  endwhile
  save, modehist,filename='/global/u1/c/clau/histsav/'+filedirectories[nfi]+'.sav'
  ioerr: close, ihistory
endfor
  
end
;**************************************************************************