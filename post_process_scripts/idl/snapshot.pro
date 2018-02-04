  pro snapshot_event,event

common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color
common snapdata,snapmenu,profile,pdf,poloidata,fluxdata,nspecies,nfield,tmax

widget_control,event.id,get_uvalue=choice
pics=0;
modedesc='drift_m=10_n=1';
snapnum=7200;  
;read,snapnum,PROMPT='snapnum='

case choice of
 
 0: begin            ; thermal ion density
	x=indgen(mpsi)
        y1=profile(*,0,0)
        y2=profile(*,1,0)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 1: begin            ; flow
	x=indgen(mpsi)
        y1=profile(*,2,0)
        y2=profile(*,3,0)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 2: begin            ; energy
	x=indgen(mpsi)
        y1=profile(*,4,0)
        y2=profile(*,5,0)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 3: begin            ; pdf in energy
	x=tmax*indgen(nvgrid)/(nvgrid-1.0)
        y1=pdf(*,0,0)
        y2=pdf(*,1,0)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 4: begin            ; pdf in pitch angle
	x=indgen(nvgrid)
        y1=pdf(*,2,0)
        y2=pdf(*,3,0)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 5: begin            ; electron density
	x=indgen(mpsi)
        y1=profile(*,0,1)
        y2=profile(*,1,1)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 6: begin            ; flow
	x=indgen(mpsi)
        y1=profile(*,2,1)
        y2=profile(*,3,1)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 7: begin            ; energy
	x=indgen(mpsi)
        y1=profile(*,4,1)
        y2=profile(*,5,1)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 8: begin            ; pdf in energy
	x=tmax*indgen(nvgrid)/(nvgrid-1.0)
        y1=pdf(*,0,1)
        y2=pdf(*,1,1)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 9: begin            ; pdf in pitch angle
	x=indgen(nvgrid)
        y1=pdf(*,2,1)
        y2=pdf(*,3,1)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 10: begin            ; EP density
	x=indgen(mpsi)
        y1=profile(*,0,2)
        y2=profile(*,1,2)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 11: begin            ; flow
	x=indgen(mpsi)
        y1=profile(*,2,2)
        y2=profile(*,3,2)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 12: begin            ; energy
	x=indgen(mpsi)
        y1=profile(*,4,2)
        y2=profile(*,5,2)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 13: begin            ; pdf in energy
	x=tmax*indgen(nvgrid)/(nvgrid-1.0)
        y1=pdf(*,0,2)
        y2=pdf(*,1,2)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 14: begin            ; pdf in pitch angle
	x=indgen(nvgrid)
        y1=pdf(*,2,2)
        y2=pdf(*,3,2)
        plot1d2p,x,y1,x,y2,"full-f","del-f"
 end

 15: begin		;phi on flux surface
	f=fluxdata(*,*,0)
	flux,f
	
	if pics eq 1 then begin
	WRITE_JPEG, strcompress('/global/u1/c/clau/snapphiflux'+string(snapnum)+modedesc+'.jpeg',/REMOVE_ALL),TVRD(),QUALITY=75
	endif
	
 end

 16: begin		;poloidal and parallel spectra                
          	f=fluxdata(*,*,0)
          	spectrum,f
	
  if pics eq 1 then begin
     WRITE_JPEG, strcompress('/global/u1/c/clau/snapspectra'+string (snapnum)+modedesc+'.jpeg',/REMOVE_ALL),TVRD(),QUALITY=75
  endif
	
 end

 17: begin		;phi on ploidal plane
	x=poloidata(*,*,nfield)
	y=poloidata(*,*,nfield+1)
	f=poloidata(*,*,0)
	poloidal,x,y,f

  ;if pics eq 1 then begin
  ;WRITE_JPEG, strcompress('/global/u1/c/clau/snappol'+string(snapnum)+modedesc+'.jpeg',/REMOVE_ALL), TVRD(),QUALITY=75
  ;endif
  
 end	

 18: begin		;radius profile of field and rms
	f=poloidata(*,*,0)
	cut1d,f,1
 end

 19: begin		;poloidal profile of field and rms
	f=poloidata(*,*,0)
	cut1d,f,2
 end

 20: begin		;a_para on flux surface
	f=fluxdata(*,*,1)
	flux,f
 end

 21: begin		;poloidal and parallel spectra                
	f=fluxdata(*,*,1)
	spectrum,f
 end

 22: begin		;a_para on ploidal plane
	x=poloidata(*,*,nfield)
	y=poloidata(*,*,nfield+1)
	f=poloidata(*,*,1)
	poloidal,x,y,f
 end	

 23: begin		;radius profile of field and rms
	f=poloidata(*,*,1)
	cut1d,f,1
 end

 24: begin		;poloidal profile of field and rms
	f=poloidata(*,*,1)
	cut1d,f,2
 end

 25: begin		;fluidne on flux surface
	f=fluxdata(*,*,2)
	flux,f
 end

 26: begin		;poloidal and parallel spectra                
	f=fluxdata(*,*,2)
	spectrum,f
 end

 27: begin		;fluidne on ploidal plane
	x=poloidata(*,*,nfield)
	y=poloidata(*,*,nfield+1)
	f=poloidata(*,*,2)
	poloidal,x,y,f
 end	

 28: begin		;radius profile of field and rms
	f=poloidata(*,*,2)
	cut1d,f,1
 end

 29: begin		;poloidal profile of field and rms
	f=poloidata(*,*,2)
	cut1d,f,2
 end

30: begin	;change poloidal & parallel mode #
	print, 'old poloidal and parallel range: m,p=',mmode,pmode
	print, 'maximal m=',mtgrid/2+1,'; maximal p=',mtoroidal/2+1
	print, 'new m,p=?'
	read, mmode,pmode
end	

31: begin	;cutoff in contour plot
	print, 'old cutoff=',cutoff
	print, 'new cutoff=?'
	read, cutoff	
end	

32: begin	;change X-window size
	size_x=1
	size_y=1
	print, 'size_x,size_y=?'
	read,size_x,size_y
        window, plotidsnap, TITLE='snapshot plot', xsize=size_x,ysize=size_y
end	

33: begin	;plot to PS
        !p.thick=2
	set_plot,'PS',/copy
	print, 'file name=?'
	title='a.ps'
	read,title
        if title ne '' then begin
	  device,file=title+'.ps',/color,Bits=8,xsize=12,ysize=12,font_size=5
	  tvlct,pseudo_color
	  truecolor=indgen(ncolor+2)
	  cvalue=truecolor(indgen(ncolor+2))
	endif
end	

34: begin            ; Quit
	set_plot,'X'
    wdelete        ; Closes plot windowse
        widget_control, event.top,/destroy
end

endcase
end
;****************************************************************

pro cut1d,f,icut

common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color

; radial profile
if icut eq 1 then begin
	x=indgen(mpsi)
	y1=fltarr(mpsi)
	y2=fltarr(mpsi)*0.0

	for i=0,mpsi-1 do begin
		y1(i)=f(0,i)

		for j=0,mtgrid-1 do begin
			y2(i)=y2(i)+f(j,i)*f(j,i)
		endfor
	endfor
	y2=sqrt(y2/mtgrid)
endif	

; poloidal profile
if icut eq 2 then begin

	x=indgen(mtgrid)
	y1=fltarr(mtgrid)
	y2=fltarr(mtgrid)*0.0

	for j=0,mtgrid-1 do begin
		y1(j)=f(j,mpsi/2)

		for i=0,mpsi-1 do begin
			y2(j)=y2(j)+f(j,i)*f(j,i)
		endfor
	endfor
	y2=sqrt(y2/mpsi)
endif

plot1d2p,x,y1,x,y2,"point value","rms"

end
;****************************************************************

pro spectrum,f
common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color

; poloidal mode
        x1=indgen(mmode)
        y1=fltarr(mmode)*0.0

        yy=fltarr(mtgrid)
        for i=0,mtoroidal-1 do begin
                yy=f(*,i)
                yy=fft(yy,1)
                y1(0)=y1(0)+(abs(yy(0)))^2
                for j=1,mmode-1 do begin
                        y1(j)=y1(j)+(abs(yy(j)))^2+(abs(yy(mtgrid-j)))^2
                endfor
        endfor
        y1=sqrt(y1/mtoroidal)
        y1=y1/mtgrid

; poloidal mode
        x2=indgen(pmode)
        y2=fltarr(pmode)*0.0

        yy=fltarr(mtoroidal)
        for i=0,mtgrid-1 do begin
                yy=f(i,*)
                yy=fft(yy,1)
                y2(0)=y2(0)+(abs(yy(0)))^2
                for j=1,pmode-1 do begin
                        y2(j)=y2(j)+(abs(yy(j)))^2+(abs(yy(mtoroidal-j)))^2
                endfor
        endfor
        y2=sqrt(y2/mtgrid)
        y2=y2/mtoroidal

	plot1d2p,x1,y1,x2,y2,"poloidal spectrum","parallel spectrum"
end
;***************************************************************

pro plot1d2p,x1,y1,x2,y2,title1,title2

common snap,plotidsnap,truecolor,ncolor,cvalue,mmode,pmode,cutoff,nvgrid,mpsi,mtgrid,mtoroidal,pseudo_color,true_color

        if !D.name eq 'X' then wset,plotidsnap
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)

        xmax=max(x1)
        xmin=min(x1)
        ymax=max(y1)
        ymin=min(y1)
        if ymin eq ymax then print, 'ymax=ymin=',ymax

        set_xy,xmin,xmax,ymin,ymax
        set_viewport,0.2,0.95,0.55,0.95
        !mtitle=title1
        plot,x1,y1,charsize=3,/xstyle,psym = -4

        xmax=max(x2)
        xmin=min(x2)
        ymax=max(y2)
        ymin=min(y2)
        if ymin eq ymax then print, 'ymax=ymin=',ymax

        set_xy,xmin,xmax,ymin,ymax
        set_viewport,0.2,0.95,0.05,0.45
        !mtitle=title2
        !noeras=1
        plot,x2,y2,charsize=3,/xstyle,psym = -4

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		truecolor=true_color
		cvalue=truecolor(indgen(ncolor+2))
                !p.thick=1
        endif
return
end
;*****************************************************************

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
;	set_viewport,0.06,0.99,0.05,0.98
        set_xy,xmin,xmax,ymin,ymax
	!x.style=1
	!y.style=1
	flevel=zmin+(zmax-zmin)*indgen(ncolor)/(ncolor-1)
  contour,f,x,y,nlevels=ncolor,c_colors=cvalue,max_value=zmax,min_value=zmin,$
		levels=flevel,/fill,charsize=2,title=''
	; plot q_min position
	!p.thick = 1
;	PLOTS, circle(1.0, 0.0, 0.173608)
 
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
pro snapshot

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
        fname = dialog_pickfile(/read, path=fpath, filter='snap*.out')
        if fname eq '' then begin
                wdelete
                return
        endif

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
        print, n_elements(profile),n_elements(pdf),n_elements(poloidata),n_elements(fluxdata)
                  
        mmode=mtgrid/5
        pmode=mtoroidal/5
	cutoff=-1.0

; make a menu
snapmenu=strarr(35)
snapmenu=["ion density","flow","energy","PDF-energy","PDF-pitch",$
     "electron density","flow","energy","PDF-energy","PDF-pitch",$
           "EP density","flow","energy","PDF-energy","PDF-pitch",$
                  "phi-flux","spectrum","poloidal","psi","theta",$
                "apara-flux","spectrum","poloidal","psi","theta",$
              "fluidne-flux","spectrum","poloidal","psi","theta",$
            "mode# range","cutoff","window size","PS file","Exit"]

xmenu,snapmenu,BASE=snapbase,SPACE=10,TITLE='snapshot panel',column=7,xpad=10,ypad=10
widget_control,snapbase,/realize
xmanager,"snapshot",snapbase,/NO_BLOCK

end
;********************************************************************************




