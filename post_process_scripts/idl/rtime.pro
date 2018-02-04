pro rtime_event,event

common rtime2d,plotidrt,truecolor,ncolor,cvalue,cutoff,nstart,nend,m0,m1,pseudo_color,true_color
common rtdata,rtmenu,mpsi,ndstep,data1di,data1de,data1df,field00,fieldrms

widget_control,event.id,get_uvalue=choice

case choice of
 
0: begin		;thermal ion particle flux
	f=data1di(*,*,0)
	rtime2d,f
 end

1: begin		;ion energy flux
	f=data1di(*,*,1)
	rtime2d,f
 end

2: begin		;electron particle flux
	f=data1de(*,*,0)
	rtime2d,f
 end

3: begin		;electron energy flux
	f=data1de(*,*,1)
	rtime2d,f
 end

4: begin		;EP particle flux
	f=data1df(*,*,0)
	rtime2d,f
 end

5: begin		;EP energy flux
	f=data1df(*,*,0)
	rtime2d,f
 end

6: begin		;zonal flow
	f=field00(*,*,0)
	rtime2d,f
 end

7: begin		;phi_rms
	f=fieldrms(*,*,0)
	rtime2d,f
 end

8: begin		;zonal current
	f=field00(*,*,1)
	rtime2d,f
 end

9: begin		;apara_rms
	f=fieldrms(*,*,1)
	rtime2d,f
 end

10: begin		;zonal fluidne
	f=field00(*,*,2)
	rtime2d,f
 end

11: begin		;fluidne_rms
	f=fieldrms(*,*,2)
	rtime2d,f
 end


12: begin	;plot ranges
	print, 'old time ranges=',nstart,nend
	print, 'old radial grid ranges=',m0,m1
	print, 'maximal time=',ndstep-1,'; maximal grid=',mpsi-1
	print, 'new time ranges=?'
	read, nstart,nend
	print, 'new grid ranges=?'
	read, m0,m1
end	

13: begin	;cutoff in contour plot
	print, 'old cutoff=',cutoff
	print, 'new cutoff=?'
	read, cutoff	
end	

14: begin	;plot to PS
        !p.thick=2
	set_plot,'PS',/copy
	print, 'file name=?'
	title='a.ps'
	read,title
        if title ne '' then begin
	      device,file=title+'.ps',/color,Bits=8,xsize=12,ysize=12,font_size=5
	      tvlct,pseudo_color
	      truecolor=indgen(ncolor+2)
	  endif
end	

15: begin            ; Quit
	set_plot,'X'
        wdelete        ; Closes plot windowse
        widget_control, event.top,/destroy
end

endcase
end
;****************************************************************
pro rtime2d,f

common rtime2d,plotidrt,truecolor,ncolor,cvalue,cutoff,nstart,nend,m0,m1,pseudo_color,true_color

        if !D.name eq 'X' then wset,plotidrt
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)
        !mtitle="data (radius,time)"

        y=indgen(m1-m0+1)
        x=indgen(nend-nstart+1)

        xmax=max(x)
        xmin=min(x)
        ymax=max(y)
        ymin=min(y)
        set_xy,xmin,xmax,ymin,ymax

; f range
	ff=f(nstart:nend,m0:m1)
        zmax=max(ff)
        zmin=min(ff)
        zmax=(0.5-cutoff)*(zmax-zmin)
        zmin=-zmax
	print,zmax,zmin

; contour plot
        set_viewport,0.1,0.9,0.1,0.9
        contour,ff,x,y,nlevels=ncolor,c_colors=cvalue,max_value=zmax,min_value=zmin,/fill,/xstyle,/ystyle,charsize=3

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		        truecolor=true_color
                !p.thick=1
        endif
return
end

;****************************************************************
pro rtime

common startup,number_plot,fpath,color_max,color_value,direct_color
common rtime2d,plotidrt,truecolor,ncolor,cvalue,cutoff,nstart,nend,m0,m1,pseudo_color,true_color
common rtdata,rtmenu,mpsi,ndstep,data1di,data1de,data1df,field00,fieldrms

plotidrt=number_plot
true_color=color_value
truecolor=true_color
pseudo_color=direct_color
ncolor=color_max
cvalue=true_color(indgen(ncolor))

; default window
  set_plot,'X'
  !p.thick=1
  !linetype=0
  !p.background=truecolor(ncolor/2)
  !p.color=truecolor(ncolor+1)
  window, plotidrt, TITLE='rtime plot', xsize=900,ysize=900

  fname="data1d.out"
  rtimeid=44
  openr, rtimeid, fname

; # of time steps
  ndstep=1
; # of radial grids: mpsi+1
  mpsi=1
; # of species: ion, electron, EP, impuries
  nspecies=1
; whether electron is loaded
  nhybrid=0
; # of variables per species: particle flux, energy flux, ...
  mpdata1d=1
; # of fields: phi, a_para, fluidne
  nfield=1
; # of variables per field: phi00, phi_rms, ...
  mfdata1d=1

  readf,rtimeid, ndstep,mpsi,nspecies,nhybrid,mpdata1d,nfield,mfdata1d
  print,ndstep,mpsi,nspecies,nhybrid,mpdata1d,nfield,mfdata1d

  ;print,'total time steps=',ndstep
  ;print,'read # of last time step'
  ;new=''
  ;read,new
  ;if new ne '' then ndstep=0+new

; data array
  data1di=fltarr(ndstep,mpsi,mpdata1d)
  data1de=fltarr(ndstep,mpsi,mpdata1d)
  data1df=fltarr(ndstep,mpsi,mpdata1d)
  field00=fltarr(ndstep,mpsi,nfield)
  fieldrms=fltarr(ndstep,mpsi,nfield)

; read rtime data
	dp=fltarr(mpsi,mpdata1d)
	df=fltarr(mpsi,nfield)

  on_ioerror,ioerr
  i=0
  while i lt ndstep do begin
;  for i=0,ndstep-1 do begin
	readf,rtimeid,dp
	data1di(i,*,*)=dp

	if nspecies eq 2 then begin	
	   if nhybrid gt 0 then begin
	      readf,rtimeid,dp
	      data1de(i,*,*)=dp
	   endif
	   if nhybrid eq 0 then begin
	      readf,rtimeid,dp
	      data1df(i,*,*)=dp
	   endif	
	endif	
	
	if nspecies eq 3 then begin
	   readf,rtimeid,dp
	   data1df(i,*,*)=dp
	   readf,rtimeid,dp
	   data1de(i,*,*)=dp
	endif	

	readf,rtimeid,df
	field00(i,*,*)=df

	readf,rtimeid,df
	fieldrms(i,*,*)=df
;  endfor		
	i=i+1
  endwhile
  
  ioerr: close, rtimeid
  ndstep=i
  nstart=0
  nend=ndstep-1
  m0=0
  m1=mpsi-1
	cutoff=-1.0

; make a menu
rtmenu=strarr(16)
rtmenu=["ion flux","energy flux",$
     "electron flux","energy flux",$
	   "EP flux","energy flux",$
                 "phi00","phi_rms",$
             "apara00","apara_rms",$
         "fluidne00","fluidne_rms",$
  "plot range","cutoff","PS file","Exit"]

xmenu,rtmenu,BASE=rtimebase,SPACE=10,TITLE='rtime panel',column=8,xpad=10,ypad=10
widget_control,rtimebase,/realize
xmanager,"rtime",rtimebase,/NO_BLOCK

end
;********************************************************************************
