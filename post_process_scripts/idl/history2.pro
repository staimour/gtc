pro history_event,event

common history,plotidh,ncolor,truecolor,nstart,nend,nfreq,ntime,xtime,tstep,pseudo_color,true_color
common hdata,hmenu,parthist,fieldhist,modehist

  widget_control,event.id,get_uvalue=choice

  case choice of
 
     0: begin		;phi	
		y=fieldhist(*,0:1,0)
		plotime,y,"phi(theta=zeta=0)","phip00(iflux_diag)"
	end

     1: begin		;RMS
		y=fieldhist(*,2:3,0)
		plotime,y,"ZF RMS","phi RMS"	
;		openw,1,'data.dat'
;		printf,1,y(0:799,0:1)
	end

     2: begin		;mode=1
		y=modehist(*,0:1,0,0)
		spectral,y,hmenu(choice)
	end

     3: begin		;mode=2
		y=modehist(*,0:1,1,0)
		spectral,y,hmenu(choice)
	end

     4: begin		;mode=3
		y=modehist(*,0:1,2,0)
		spectral,y,hmenu(choice)
	end

     5: begin		;mode=4
		y=modehist(*,0:1,3,0)
		spectral,y,hmenu(choice)
	end

     6: begin		;mode=5
		y=modehist(*,0:1,4,0)
		spectral,y,hmenu(choice)
	end

     7: begin		;mode=6
		y=modehist(*,0:1,5,0)
		spectral,y,hmenu(choice)
	end

     8: begin		;mode=7
		y=modehist(*,0:1,6,0)
		spectral,y,hmenu(choice)
	end

     9: begin		;mode=8
	
		y=modehist(*,0:1,7,0)
		spectral,y,hmenu(choice)
	end

     10: begin		;a_par	
		y=fieldhist(*,0:1,1)
		plotime,y,"a_par(theta=zeta=0)","a_par00(iflux_diag)"
	end

     11: begin		;zonal fluidne
		y=fieldhist(*,2:3,1)
		plotime,y,"ZF RMS","a_par RMS"
	end

     12: begin		;mode=1
		y=modehist(*,0:1,0,1)
		spectral,y,hmenu(choice)
	end

     13: begin		;mode=2
		y=modehist(*,0:1,1,1)
		spectral,y,hmenu(choice)
	end

     14: begin		;mode=3
		y=modehist(*,0:1,2,1)
		spectral,y,hmenu(choice)
	end

     15: begin		;mode=4
		y=modehist(*,0:1,3,1)
		spectral,y,hmenu(choice)
	end

     16: begin		;mode=5
		y=modehist(*,0:1,4,1)
		spectral,y,hmenu(choice)
	end

     17: begin		;mode=6
		y=modehist(*,0:1,5,1)
		spectral,y,hmenu(choice)
	end

     18: begin		;mode=7
		y=modehist(*,0:1,6,1)
		spectral,y,hmenu(choice)
	end

     19: begin		;mode=8
	
		y=modehist(*,0:1,7,1)
		spectral,y,hmenu(choice)
	end

     20: begin		;fluidne	
		y=fieldhist(*,0:1,2)
		plotime,y,"fluidne(theta=zeta=0)","fluidne_00(iflux_diag)"
	end

     21: begin		;zonal fluidne
		y=fieldhist(*,2:3,2)
		plotime,y,"ZF RMS","fluidne RMS"
	end

     22: begin		;mode=1
		y=modehist(*,0:1,0,2)
		spectral,y,hmenu(choice)
	end

     23: begin		;mode=2
		y=modehist(*,0:1,1,2)
		spectral,y,hmenu(choice)
	end

     24: begin		;mode=3
		y=modehist(*,0:1,2,2)
		spectral,y,hmenu(choice)
	end

     25: begin		;mode=4
		y=modehist(*,0:1,3,2)
		spectral,y,hmenu(choice)
	end

     26: begin		;mode=5
		y=modehist(*,0:1,4,2)
		spectral,y,hmenu(choice)
	end

     27: begin		;mode=6
		y=modehist(*,0:1,5,2)
		spectral,y,hmenu(choice)
	end

     28: begin		;mode=7
		y=modehist(*,0:1,6,2)
		spectral,y,hmenu(choice)
	end

     29: begin		;mode=8
	
		y=modehist(*,0:1,7,2)
		spectral,y,hmenu(choice)
	end

     30: begin		;ion density & entropy	
		y=parthist(*,0:1,0)
		plotime,y,"density delta_f","entropy delta_f^2"	
	end

     31: begin		;momentum
		y=parthist(*,2:3,0)
		plotime,y,"paralle flow u", "delta_u"
	end

     32: begin		;energy
		y=parthist(*,4:5,0)
		plotime,y,"energy E-1.5","delta_E"
	end

     33: begin		;particle & momentum flux
		y=parthist(*,6:7,0)
		plotime,y,"particle flux","momentum flux"
	end

     34: begin		;energy flux
		y=parthist(*,8:9,0)
		plotime,y,"energy flux","total density"
	end

     35: begin		;EP density & entropy	
		y=parthist(*,0:1,2)
		plotime,y,"density delta_f","entropy delta_f^2"	
	end

     36: begin		;momentum
		y=parthist(*,2:3,2)
		plotime,y,"paralle flow u", "delta_u"	
	end

     37: begin		;energy
		y=parthist(*,4:5,2)
		plotime,y,"energy E-1.5","delta_E"
	end

     38: begin		;particle & momentum flux
		y=parthist(*,6:7,2)
		plotime,y,"particle flux","momentum flux"
	end

     39: begin		;energy flux
		y=parthist(*,8:9,2)
		plotime,y,"energy flux","total density"
	end

     40: begin		;electron density & entropy	
		y=parthist(*,0:1,1)
		plotime,y,"density delta_f","entropy delta_f^2"		
	end

     41: begin		;momentum
		y=parthist(*,2:3,1)
		plotime,y,"paralle flow u", "delta_u"	
	end

     42: begin		;energy
		y=parthist(*,4:5,1)
		plotime,y,"energy E-1.5","delta_E"
	end

     43: begin		;particle & momentum flux
		y=parthist(*,6:7,1)
		plotime,y,"particle flux","momentum flux"
	end

     44: begin		;energy flux
		y=parthist(*,8:9,1)
		plotime,y,"energy flux","total density"
	end

	45: begin       ;plot range
		print, 'current range=',nstart,',',nend,'.   maximal nend=0,',ntime-1
        	print, 'new range: nstart,nend=?'
		oldend=nend
        	read, nstart,nend
		if nend gt oldend then nend=oldend
	end

	46: begin       ;change frequency range
		print, 'current nfreq=',nfreq
        	print, 'new nfreq=?'
        	read, nfreq
	end

	47: begin       ;change window size
        	size_x=1
        	size_y=1
        	print, 'size_x,size_y=?'
        	read,size_x,size_y
        	window, plotidh, TITLE='history', xsize=size_x,ysize=size_y
	end

	48: begin       ;plot to PS
		!p.thick=2
        	print, 'file name=?'
        	title='a'
        	read,title
        	set_plot,'PS'
        	device,file=title+'.ps',/color,Bits_per_pixel=8,xsize=16,ysize=12,font_size=5
		tvlct,pseudo_color
		truecolor=indgen(ncolor+2)
	end
		
	49: begin            ; Quit
        	wdelete        ; Closes plot windows
         	widget_control, event.top,/destroy
 	end

endcase
end

;******************************************************************
pro plotime,y,title0,title1

common history,plotidh,ncolor,truecolor,nstart,nend,nfreq,ntime,xtime,tstep,pseudo_color,true_color

	if !D.name eq 'X' then wset,plotidh
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)

	xx=xtime(nstart:nend)
	y0=y(nstart:nend,0)
	y1=y(nstart:nend,1)

	xmax=max(xx)
        xmin=min(xx)
	ymax=max(y0)
        ymin=min(y0)
	if ymin eq ymax then print, 'ymax=ymin=',ymax

        set_xy,xmin,xmax,ymin,ymax
	set_viewport,0.2,0.95,0.55,0.95
	!mtitle=title0
	plot,xx,y0,charsize=3

	ymax=max(y1)
        ymin=min(y1)
	if ymin eq ymax then print, 'ymax=ymin=',ymax

        set_xy,xmin,xmax,ymin,ymax
	set_viewport,0.2,0.95,0.05,0.45
	!mtitle=title1
        !noeras=1
	plot,xx,y1,charsize=3

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		truecolor=true_color
		!p.thick=1
        endif
return
end		

;******************************************************************
pro spectral,y,title

common history,plotidh,ncolor,truecolor,nstart,nend,nfreq,ntime,xtime,tstep,pseudo_color,true_color
;common idlh,plotidh,ncolor,np,nm,glength,truecolor

        if !D.name eq 'X' then wset,plotidh
	!mtitle=title
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)

; panel 1: mode history of real and imaginary components
	yr=y(nstart:nend,0)
	yi=y(nstart:nend,1)
	xx=xtime(nstart:nend)

	xmax=max(xx)
        xmin=min(xx)
        ymax=max([yr,yi])
        ymin=min([yr,yi])
        set_xy,xmin,xmax,ymin,ymax
        set_viewport,0.14,0.49,0.6,0.95
        !linetype=0
        !p.color=truecolor(ncolor+1)
        plot,xx,yr,charsize=2,/ystyle
        !p.color=truecolor(ncolor*3/4)
        oplot,xx,yi

; panel 2: mode amplitude history
	ya=sqrt(yr*yr+yi*yi)
	ypow=alog10(ya)

        ymax=max(ypow)
        ymin=min(ypow)
        set_xy,xmin,xmax,ymin,ymax
        set_viewport,0.63,0.98,0.6,0.95
        !noeras=1
        !p.color=truecolor(ncolor+1)
        !mtitle="mode amplitude"
        plot,xx,ypow,charsize=2,/ystyle

; panel 3: mode amplitude normalized by growth rate
        gamma=alog(ya(nend-nstart)/ya(0))/(nend-nstart)
        xpow=indgen(nend-nstart)
	yr=yr/exp(gamma*xpow)
	yi=yi/exp(gamma*xpow)
        print,"growth rate=",gamma/tstep
                
        ymin=min([yr,yi])
        ymax=max([yr,yi])
        set_xy,xmin,xmax,ymin,ymax
        set_viewport,0.14,0.49,0.1,0.45
        !noeras=1
	!mtitle="normilized by gamma"
	!p.color=truecolor(ncolor+1)
        plot,xx,yr,charsize=2,/ystyle
        !p.color=truecolor(ncolor*3/4)
        oplot,xx,yi

; panel 4: power spectral: phi=exp[-i(omega*t+m*theta-n*zeta)]
	power=complex(yr,yi)
        power=fft(power,-1)
        ypow=abs(power)
        yp=fltarr(2*nfreq)
        xp=fltarr(2*nfreq)
        for i=0,nfreq-2 do begin
            yp(i)=ypow(i+nend-nstart-nfreq+1)
            xp(i)=(i-nfreq+1)*6.283185/((nend-nstart)*tstep)
        end
        for i=0, nfreq do begin
            yp(nfreq-1+i)=ypow(i)
            xp(nfreq-1+i)=i*6.283185/((nend-nstart)*tstep)
        end

        xmax=max(xp)
        xmin=min(xp)
        ymax=max(yp)
        ymin=min(yp)
        set_xy,xmin,xmax,ymin,ymax
        set_viewport,0.63,0.98,0.1,0.45
        !mtitle="frequency"
        !noeras=1
        !p.color=truecolor(ncolor+1)
        plot,xp,yp,/xstyle,charsize=2,/ystyle

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		truecolor=true_color
		!p.thick=1
        endif
return
end
;****************************************************************
pro history
	
common startup,number_plot,fpath,color_max,color_value,direct_color
common history,plotidh,ncolor,truecolor,nstart,nend,nfreq,ntime,xtime,tstep,pseudo_color,true_color
common hdata,hmenu,parthist,fieldhist,modehist

plotidh=number_plot
ncolor=color_max
true_color=color_value
truecolor=color_value
pseudo_color=direct_color

; default window
  	set_plot,'X'
	!p.thick=1
        !linetype=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)
        window, plotidh, TITLE='history-plot', xsize=900,ysize=900

	fname="history.out"

	ihistory=11
	openr, ihistory, fname

; # of time steps
	ndstep=1
; # of species: ion, electron, EP, impuries
	nspecies=1
; # of quantities per species: density,entropy,momentum,energy, fluxes
	mpdiag=1	
; # of field variables: phi, a_par, fluidne
	nfield=1
; # of modes per field: (n,m)
	modes=1
; # of quantities per field: rms, single spatial point, zonal components
	mfdiag=1
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
	i=0
	while i lt ndstep do begin
;	for i=0,ndstep-1 do begin
		readf,ihistory,partdata,fieldtime,fieldmode
		parthist(i,*,*)=partdata
		fieldhist(i,*,*)=fieldtime
		modehist(i,*,*,*)=fieldmode
;	print,i+1,modehist(i,0,4,2),modehist(i,1,4,2)
;	endfor
		i=i+1
	endwhile
	ioerr: close, ihistory
	ndstep=i
	
	ntime=ndstep
	xtime=indgen(ntime)
	nstart=0
	nend=ntime-1
	nfreq=(nend-nstart)/10

;	np=(ntime-nstart)/50
;	if np lt 2 then np=2
;	nm=(ntime-nstart)/4

hmenu=strarr(50)
hmenu=["phi","RMS","mode1","mode2","mode3","mode4","mode5","mode6","mode7","mode8",$
     "apara","RMS","mode1","mode2","mode3","mode4","mode5","mode6","mode7","mode8",$
   "fluidne","RMS","mode1","mode2","mode3","mode4","mode5","mode6","mode7","mode8",$
         "ion density","momentum","energy","pmflux","eflux",$	
          "EP density","momentum","energy","pmflux","eflux",$	
    "electron density","momentum","energy","pmflux","eflux",$	
      "time range","frequency range","window size","PS file","Exit"]

xmenu,hmenu,BASE=hbase,SPACE=5,TITLE='history-panel',column=5,xpad=8,ypad=5
widget_control,hbase,/realize
xmanager,"history",hbase,/NO_BLOCK

end
;**************************************************************************
