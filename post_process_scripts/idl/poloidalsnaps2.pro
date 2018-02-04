pro poloidal,x,y,f

  common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color
  
  if !D.name eq 'X' then wset,plotidsnap
  !noeras=0
  !p.background=truecolor(ncolor/2)
  !p.color=truecolor(ncolor+1)
  
  ; 8-bits color
  ;        !p.background=ncolor
  ;        !p.color=ncolor+1
  
  ; ploting range
  length=1.0
  x=x*length
  y=y*length
  xmax=max(x)
  xmin=min(x)
  ymax=max(y)
  ymin=min(y)
  ;zmax=max(f)
  zmax=max(abs(f))
  ;zmin=min(f)
  zmin=-zmax
  if zmax eq zmin then print,zmax
  
  ; contour plot
  xw=0.1+0.9*(xmax-xmin)/(ymax-ymin)
  set_viewport,0.1,xw,0.1,0.9
  ; set_viewport,0.06,0.99,0.05,0.98
  set_xy,xmin,xmax,ymin,ymax
  !x.style=1
  !y.style=1
  flevel=zmin+(zmax-zmin)*indgen(ncolor)/(ncolor-1)
  contour,f,x,y,nlevels=ncolor,c_colors=cvalue,max_value=zmax,min_value=zmin,$
    levels=flevel,/fill,charsize=2,title=''
  ; plot q_min position
  !p.thick = 1
  ; PLOTS, circle(1.0, 0.0, 0.173608)
  
  if !D.name eq 'PS' then begin
    device,/close
    set_plot,'X'
    truecolor=true_color
    cvalue=truecolor(indgen(ncolor+2))
    !p.thick=1
  endif
  return
end
;****************************************************************

pro flux,f

  common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color
  
  if !D.name eq 'X' then wset,plotidsnap
  !noeras=0
  !p.background=truecolor(ncolor/2)
  !p.color=truecolor(ncolor+1)
  !mtitle="flux surface (parallel,poloidal)"
  
  x=2.0*3.1415926*((indgen(mtgrid)+1.0)/mtgrid)
  y=2.0*3.1415926*((indgen(mtoroidal)+1.0)/mtoroidal)
  ff=fltarr(mtoroidal,mtgrid)
  for i=0,mtoroidal-1 do begin
    ff(i,*)=f(*,i)
  endfor
  
  xmax=max(x)
  xmin=min(x)
  ymax=max(y)
  ymin=min(y)
  ;        del=(ymax-ymin-xmax+xmin)/2.0
  ;        if del gt 0.0 then xmax=xmax+0.8*del
  ;        if del gt 0.0 then xmin=xmin-0.8*del
  ;        if del lt 0.0 then ymax=ymax-0.8*del
  ;        if del lt 0.0 then ymin=ymin+0.8*del
  set_xy,xmin,xmax,ymin,ymax
  
  ; ff range
  zmax=max(ff)
  zmin=min(ff)
  zmax=(0.5-cutoff)*(zmax-zmin)
  zmin=-zmax
  
  ; contour plot
  set_viewport,0.1,0.9,0.1,0.9
  contour,ff,y,x,nlevels=ncolor,c_colors=cvalue,max_value=zmax,min_value=zmin,/fill,/xstyle,/ystyle,charsize=3
  
  if !D.name eq 'PS' then begin
    device,/close
    set_plot,'X'
    truecolor=true_color
    cvalue=truecolor(indgen(ncolor+2))
    !p.thick=1
  endif
  return
end
;****************************************************************
pro poloidalsnaps

  common startup,number_plot,fpath,color_max,color_value,direct_color
  common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color
  common snapdata,snapmenu,profile,pdf,poloidata,fluxdata,nspecies,nfield,tmax
  
  plotidsnap=number_plot
  true_color=color_value
  truecolor=true_color
  pseudo_color=direct_color
  ncolor=color_max
  cvalue=truecolor(indgen(ncolor+2))
  
  ; default window
  set_plot,'X'
  !p.thick=1
  !linetype=0
  !p.background=truecolor(ncolor/2)
  !p.color=truecolor(ncolor+1)
  window, plotidsnap, TITLE='snapshot plot', xsize=920,ysize=920
  fname="snap00000.out"
  startat=0
  snapsize=250
  numofsnaps=55
for i=0,numofsnaps do begin
  
  snapnumname=strcompress(string(i*snapsize+startat),/REMOVE_ALL)
  add=5-strlen(snapnumname);
  if add gt 0 then snapnumname=string(bytarr(add)+48B)+snapnumname
  if i eq 0 then snapnumname='00001'  
  fname=strcompress("snap"+snapnumname+".out",/REMOVE_ALL)
  
  snapid=22
  openr, snapid, fname
  
  ; # of species: ion, electron, EP, impuries
  nspecies=1
  ; # of field variables: phi, a_para, fluidne
  nfield=1
  ; # of grids in energy and pitch
  nvgrid=1
  ; # of radial grids: mpsi+1
  mpsi=1
  ; # of poloidal grids
  mtgrid=1
  ; # of toroidal (parallel) grids
  mtoroidal=1
  ; upper bound of temperature
  tmax=1.0
  
  readf,snapid, nspecies,nfield,nvgrid,mpsi,mtgrid,mtoroidal,tmax
  print, nspecies,nfield,nvgrid,mpsi,mtgrid,mtoroidal,tmax
  
  ; data array
  profile=fltarr(mpsi,6,nspecies)
  pdf=fltarr(nvgrid,4,nspecies)
  poloidata=fltarr(mtgrid,mpsi,nfield+2)
  fluxdata=fltarr(mtgrid,mtoroidal,nfield)
  ; read snapshot data
  readf,snapid,profile,pdf,poloidata,fluxdata
  close, snapid
  
  mmode=mtgrid/5
  pmode=mtoroidal/5
  cutoff=-1.0
  
  ; make a menu
  
  x=poloidata(*,*,nfield)
  y=poloidata(*,*,nfield+1)
  f=poloidata(*,*,0)
  poloidal,x,y,f
  
  image3d = TVRD(TRUE=1)
  WRITE_JPEG, strcompress('/global/u1/c/clau/edge_turnsoff/'+snapnumname+'.jpeg',/REMOVE_ALL), image3d, TRUE=1, QUALITY=75
  print, i, snapnumname
  wait,1.0

 endfor
  
end
;********************************************************************************




