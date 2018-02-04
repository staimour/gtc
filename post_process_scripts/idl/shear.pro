pro shear
common startup,number_plot,fpath,color_max,color_value,direct_color

  ncolor=color_max
  truecolor=color_value

; read shear and flux data
  nterm=1
  mpsi=1
  openr, 1, 'sheareb.out'
  readf,1,nterm
  readf,1,mpsi

  print,'# of radial grid=',mpsi,'  # of ntime step=?'
  ntimeall=1
  read,ntimeall
  sdata=fltarr(mpsi,nterm,ntimeall)
  readf,1,sdata
  close,1

  nbegin=0
  ntime=ntimeall
  y=indgen(mpsi)
  time=nbegin+indgen(ntime)
  f1=fltarr(ntime,mpsi)
  f2=fltarr(ntime,mpsi)
  for j=0,ntime-1 do begin
      for i=0,mpsi-1 do begin
          f1(j,i)=sdata(i,0,j+nbegin)*sdata(i,0,j+nbegin)
          f2(j,i)=sdata(i,1,j+nbegin)
      end
  end

  set_plot,'X'
  layer=ncolor
  cvalue=color_value(indgen(layer)*ncolor/layer)
  window, 1, xsize=1600,ysize=800  
  
  set_viewport,0.1,0.5,0.1,0.9
  zmax=max(f1)
  zmin=min(f1)
  flevel=zmin+(zmax-zmin)*indgen(layer)/(layer-1)
  contour,f1,time,y,nlevels=layer,c_colors=cvalue,max_value=zmax,min_value=zmin,levels=flevel,/fill

  read,dump

  !noeras=1
  set_viewport,0.55,0.95,0.1,0.9
  zmax=max(f2)
  zmin=min(f2)
  flevel=zmin+(zmax-zmin)*indgen(layer)/(layer-1)
  contour,f2,time,y,nlevels=layer,c_colors=cvalue,max_value=zmax,min_value=zmin,levels=flevel,/fill

end
