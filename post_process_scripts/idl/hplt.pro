pro history_event,event

common idlh,plotidh,ncolor,np,nm,glength,truecolor
common data,f,hname,xtime,nbegin,naverage,ntime,delt

  widget_control,event.id,get_uvalue=choice

  case choice of 
    0: begin            ; Quit
         wdelete        ; Closes plot windows
         widget_control, event.top,/destroy
       end

    1: begin		;empty

; subtract noise for WTT
;		y=f(13,*)*glength*glength
		y=f(13,*)*glength*glength-f(15,*)*20.0
		!mtitle='scaling'
		yy=f(22,*)*glength
		
		y=0.95*y/max(y)
		yy=0.95*yy/max(yy)

		y5=fltarr(5,ntime)
		y5(0,*)=y
		y5(1,*)=yy
		y5(2,*)=yy/y*0.2
		y5(3,*)=sqrt(f(12,*))*glength
		ymax=max(y5(3,*))
		y5(3,*)=0.95*y5(3,*)/ymax
		y5(4,*)=0.0
		plot5,xtime,y5,nbegin,ntime	
;	xx=xtime(0:3249)
;	yy=y5(0:4,0:3249)
;	for i=1,750 do begin
;		yy(0:4,2499+i)=y5(0:4,2499+2*i)
;	end
;		plot5,xx,yy,0,3249
	end

    2: begin		; orbit_poloidal
		x=f(0,*)
		y=f(1,*)
                
;                openw,2,'phirms.dat'
;                printf,2,ntime
;                printf,2,x,y
;                close,2 
                print,x
		!mtitle=hname(choice)
		plot1,x,y,nbegin,ntime
       end

     3: begin		;orbit_flux
		y=f(2,*)
		x=f(3,*)
	
		!mtitle=hname(choice)
		plot1,x,y,nbegin,ntime		
	end

     4: begin		;weight
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end


     5: begin		;u_para
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end
	
     6: begin		;energy
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     7: begin		;angular_momentum
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     8: begin		;ddeni
		y=f(choice,*)
;		y1=y(nbegin:ntime-1)
;		y2=y(nbegin:ntime-1)
;		for i=1,naverage do begin
;			y1=shift(y(nbegin:ntime-1),1)
;			y1(ntime-1-nbegin)=y(ntime-nbegin)
;			y2=shift(y(nbegin:ntime-1),-1)
;			y2(0)=y(nbegin)
;			y(nbegin:ntime-1)=0.5*y(nbegin:ntime-1)+0.25*(y1+y2)
;		endfor
		
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     9: begin		;ddene
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime
	end

     10: begin		;radial mode
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     11: begin		;field energy
	
		y=f(choice,*)
                y=y
;		y=y/y(0)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime

;	openw,1,'flow.dat'
;	printf,1,y
;	close,1

;	xx=xtime(nbegin:ntime-1)
;	yy=f(choice-1,nbegin:ntime-1)
;	oplot,xx,yy

	end

     12: begin		;entropyi
; data: RMS of V_zonal/V_th
; normalize zonal flow by diamagnetic velocity	
		y=sqrt(f(choice,*))
		y=f(choice,*)
		!mtitle=hname(choice)
		xlog=alog(xtime)
		ylog=alog10(y)
;		plot1,xtime,y*glength*glength,nbegin,ntime		
		plot1,xtime,y,nbegin,ntime		
;	growth
;	dy=y*0.0
;	nd=1
;	for n=nd,ntime-nd-1 do begin
;		dy(n)=(y(n+nd)-y(n-nd))/y(n)
;	endfor
;	dy=0.5*dy/(2.0*nd*delt)
; average
;	ady=0.0
;	for n=nbegin-1, ntime-1 do begin
;	        ady=ady+dy(n)
;	endfor	
;	ady=ady/(ntime-nbegin)
;        print,'growth rate=',ady

;	plot1,xtime,dy,nbegin,ntime	

	end

     13: begin		;entropye
; subtract noise for WTT
;		y=f(choice,*)*glength*glength-f(15,*)*20.0
; data: RMS of e*phi/T
; mixing length nornalization: 1/rho*
		y=f(choice,*)
		!mtitle=hname(choice)
;		plot1,xtime,alog10(y),nbegin,ntime	
		plot1,xtime,y,nbegin,ntime	
;                gamma=alog(y(ntime-1)/y(nbegin-1))/((ntime-nbegin)*delt)
;                print,"growth rate=",gamma
	end

     14: begin		;dlfowi
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     15: begin		;dflowe
	
		y=f(choice,*)

; Monte-Carlo calculation
;	y=y+f(choice-3,*)*f(choice-3,*)+f(choice-2,*)*f(choice-2,*)
;	y=y+f(choice-3,*)*f(choice-3,*)
;	diff=y
;	for i=nbegin+1,ntime-2 do begin 
;		diff(i)=102.2*(y(i+1)-y(i-1))/(4.0*delt)
;		diff(i)=97.68*(y(i+1)-y(i-1))/(4.0*delt)
;		diff(i)=105.0*(y(i+1)-y(i-1))/(4.0*delt)
;		diff(i)=130.2*(y(i+1)-y(i-1))/(4.0*delt)
;		diff(i)=138.0*(y(i+1)-y(i-1))/(4.0*delt)
;		diff(i)=172.0*(y(i+1)-y(i-1))/(4.0*delt)
;	endfor

;	diff(nbegin)=diff(nbegin+1)
;	diff(ntime-1)=diff(ntime-2)
	
;	diff=1.0-f(choice-7,*)/diff
;	diff=diff-f(choice-7,*)

		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
;		plot1,xtime,diff,nbegin,ntime	
	end

     16: begin		;pfluxi
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     17: begin		;pfluxe
	
		y=f(choice,*)
		!mtitle=hname(choice)
		plot1,xtime,y,nbegin,ntime	
	end

     18: begin		;efluxi
	
		y=f(choice,*)
		y1=y(nbegin:ntime-1)
		y2=y(nbegin:ntime-1)
;		for i=1,naverage do begin
;			y1=shift(y(nbegin:ntime-1),1)
;			y1(ntime-1-nbegin)=y(ntime-nbegin)
;			y2=shift(y(nbegin:ntime-1),-1)
;			y2(0)=y(nbegin)
;			y(nbegin:ntime-1)=0.5*y(nbegin:ntime-1)+0.25*(y1+y2)
;		endfor
		!mtitle=hname(choice)
;                y=y*317.0
		plot1,xtime,y,nbegin,ntime	
;openw,11,'chi-i.dat'
;printf,11,y
;close,11
	end

     19: begin		;efluxe
	
		!mtitle=hname(choice)
		y=f(choice,*)
 ;               y=y*317.0
		plot1,xtime,y,nbegin,ntime
                
;openw,11,'chi-e.dat'
;printf,11,y
;close,11

; Monte-Carlo calculation
;	diff=y
;	for i=nbegin+1,ntime-2 do begin 
;		diff(i)=1.5*(y(i+1)-y(i-1))/(4.0*delt)/2.0
;	endfor

;	diff(nbegin)=diff(nbegin+1)
;	diff(ntime-1)=diff(ntime-2)
;	diff=diff/f(choice+3,*)
;		plot1,xtime,diff,nbegin,ntime	
	end

     20: begin		;radial profile of energy flux
         
		y5=f(choice:choice+4,*)*glength
		!mtitle=hname(choice)
		plot5,xtime,y5,nbegin,ntime		
	end

     21: begin		;mode=1
	
		ydum=f(choice+20,*)
		yy=f(choice+4:choice+5,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt
	end

     22: begin		;mode=2

		ydum=f(choice+21,*)
		yy=f(choice+5:choice+6,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt
	end

     23: begin		;mode=3
	
		ydum=f(choice+22,*)
		yy=f(choice+6:choice+7,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt
	end

     24: begin		;mode=4

		ydum=f(choice+23,*)
		yy=f(choice+7:choice+8,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt
	end

     25: begin		;mode=5

		ydum=f(choice+24,*)
		yy=f(choice+8:choice+9,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt
	end
     26: begin		;mode=6

		yy=f(choice+9:choice+10,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
     27: begin		;mode=7
	
		yy=f(choice+10:choice+11,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
     28: begin		;mode=8
	
		yy=f(choice+11:choice+12,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
	
     29: begin		;mode=1
	
		yy=f(choice+12:choice+13,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end

     30: begin		;mode=2

		yy=f(choice+13:choice+14,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end

     31: begin		;mode=3
	
		yy=f(choice+14:choice+15,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end

     32: begin		;mode=4

		yy=f(choice+15:choice+16,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt
	end

     33: begin		;mode=5

		yy=f(choice+16:choice+17,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
     34: begin		;mode=6

		yy=f(choice+17:choice+18,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
     35: begin		;mode=7
	
		yy=f(choice+18:choice+19,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
     36: begin		;mode=8
	
		yy=f(choice+19:choice+20,*)
		!mtitle=hname(choice)
		spectrum,xtime,yy,nbegin,ntime,delt	
	end
37: begin		;peak
	nbin='3'
	print,'which radial bin to plot?'
	read,nbin
	if nbin	eq '' then nbin='3'
	i=19+nbin
	y=f(i,*)
;	y=(f(i,*)+f(i-1,*)+f(i+1,*))/3
;	y=y*glength
	y=y*125.0
;*glength
;	y=(f(i,*)+f(i-1,*)+f(i+1,*))*glength/3.0
;	y=f(i,*)/f(13,*)^2*glength
;	average=0.0
;	for j=750,1149 do begin
;	for j=1250,1749 do begin
;		average=average+y(j)
;	endfor
;	print,average/500.0

	!mtitle=hname(choice)
	plot1,xtime,y,nbegin,ntime

;        openw,1,'chi.dat'
;        printf,1,y
;        close,1
end

38: begin       ;starting time
	print, 'current nbegin, naverage=',nbegin, naverage
        print, 'nbegin, naverage=?'
        read, nbegin,naverage
	if nm gt ntime-nbegin then nm=ntime-nbegin
end

39: begin       ;change frequency range
	print, 'current np, nm=',np, nm
        print, 'new np,nm=?'
        read, np,nm
end

40: begin       ;plot to X window
        window, plotidh, TITLE='history', xsize=1000,ysize=500
end

41: begin       ;plot to X window
        size_x=1
        size_y=1
        print, 'size_x,size_y=?'
        read,size_x,size_y
        window, plotidh, TITLE='history', xsize=size_x,ysize=size_y

end

42: begin       ;plot to PS

	!p.thick=2
        set_plot,'PS'
        print, 'file name=?'
        title='a'
        read,title
        device,file=title+'.ps',/color,Bits_per_pixel=8,xsize=16,ysize=12,font_size=5
end
	
  endcase
  end

;******************************************************************
; single line plot
pro plot1,x,y,nbegin,ntime
common idlh,plotidh,ncolor,np,nm,glength,truecolor,pseudo_color

	if !D.name eq 'X' then wset,plotidh
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)
        !p.color=truecolor(ncolor+1)

	set_viewport,0.15,0.95,0.1,0.9

;	xx=x(nbegin:ntime-1)-x(nbegin)
	xx=x(nbegin:ntime-1)
	yy=y(nbegin:ntime-1)
;	xx=(x(nbegin:ntime-1)+1.0)*4000.0/0.358
;	yy=y(nbegin:ntime-1)*4000.0/0.358

;	average=total(yy)/(ntime-nbegin)
;	print,'average=',average
;	for i=nbegin, ntime-1 do begin
;		print, i,y(i)
;	endfor

	xx=xx+1.0
	xmax=max(xx)
        xmin=min(xx)
	ymax=max(yy)
        ymin=min(yy)
;	print,ymax
;	ymax=22
;	ymin=-ymax
;	xmax=2500
;	ymin=0.0
	if ymin eq ymax then print, 'y=',ymax

        set_xy,xmin,xmax,ymin,ymax
;	plot,xx,yy,xrange=[xmin,xmax],yrange=[ymin,ymax]
;	plot,charsize=1.5,[0,xx],[0,yy],/ystyle
	plot,charsize=1.5,xx,yy,/ystyle

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		truecolor=true_color
		!p.thick=1
        endif

;                openw,1,'data.dat'
;                printf,1,yy
;                close,1

return
end		
;******************************************************************
; plot 5 lines
pro plot5,xa,ya5,nbegin,ntime
common idlh,plotidh,ncolor,np,nm,glength,truecolor,pseudo_color

	if !D.name eq 'X' then wset,plotidh
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)
	set_viewport,0.15,0.95,0.1,0.9

	x=xa(nbegin:ntime-1)-xa(nbegin)
	y5=ya5(*,nbegin:ntime-1)

	xmax=max(x)
        xmin=min(x)
	ymax=max(y5)
        ymin=min(y5)
	print,y5(*,ntime-1-nbegin)	
	
;	xmax=150.0	
;	ymin=0.0
	
        set_xy,xmin,xmax,ymin,ymax

	y=y5(0,*)
	plot,x,y
	xyouts,0.2,0.9,charsize=2.0,'r=1',/normal

        !p.color=truecolor(ncolor*3/8)
	y=y5(1,*)
	oplot,x,y
	xyouts,0.2,0.8,charsize=2.0,'r=2',/normal

        !p.color=truecolor(ncolor*13/16)
	y=y5(2,*)
	oplot,x,y
	xyouts,0.2,0.7,charsize=2.0,'r=3',/normal

        !p.color=truecolor(ncolor/8)
	y=y5(3,*)
	oplot,x,y
	xyouts,0.2,0.6,charsize=2.0,'r=4',/normal

        !p.color=truecolor(ncolor*5/8)
	y=y5(4,*)
	oplot,x,y
	xyouts,0.2,0.5,charsize=2.0,'r=5',/normal

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		truecolor=true_color
		!p.thick=1
        endif

return
end		
;******************************************************************
; single spectrum
pro spectrum,x,yy,nbegin,ntime,delt
common idlh,plotidh,ncolor,np,nm,glength,truecolor,pseudo_color

        if !D.name eq 'X' then wset,plotidh
        !noeras=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)
        set_viewport,0.1,0.5,0.55,0.95

; panel 1: mode history of real and imaginary components
		yr=fltarr(ntime-nbegin)
		yi=fltarr(ntime-nbegin)
		ya=fltarr(ntime-nbegin)
		for i=0,ntime-nbegin-1 do begin
			yr(i)=yy(0,i+nbegin)
			yi(i)=yy(1,i+nbegin)

		endfor
;		yi(0:599)=yi(0:599)/0.1483640401E-02
		ya=sqrt(yr*yr+yi*yi)
                xpow=x(nbegin:ntime-1)-x(nbegin)
                xmax=max(xpow)
                xmin=min(xpow)
;		xmin=0.0
                ymax=max([yr,yi])
                ymin=min([yr,yi])
                set_xy,xmin,xmax,ymin,ymax
                !linetype=0
                plot,xpow,yr
;                xyouts,0.1,0.9,charsize=2.0,'real',/normal

                !p.color=truecolor(ncolor*3/4)
                oplot,xpow,yi
;		print,yi
;                xyouts,0.1,0.86,charsize=2.0,'imaginary',/normal

; panel 2: mode amplitude history
        !noeras=1
;        set_viewport,0.1,0.9,0.1,0.9
        set_viewport,0.55,0.95,0.55,0.95

                !p.color=truecolor(ncolor+1)
                !mtitle="mode amplitude .vs. t"
		ypow=alog10(ya)
;		ypow=ya
                ymax=max(ypow)
                ymin=min(ypow)
;                ymax=0.01
                set_xy,xmin,xmax,ymin,ymax
                plot,xpow,ypow

; panel 3: mode amplitude normalized by growth rate
        !noeras=1
        set_viewport,0.14,0.54,0.05,0.45

                !mtitle="mode"
                gamma=alog(ya(ntime-nbegin-1)/ya(0))/((ntime-nbegin)*delt)
                print,"growth rate=",gamma
;                gamma=0.0
                xpow=indgen(ntime-nbegin)*delt
		yr=yr/exp(gamma*xpow)
		yi=yi/exp(gamma*xpow)
                xpow=x(nbegin:ntime-1)-x(nbegin)
; subtract average
		ysize=size(yr)
		mean=total(yr)/ysize(1)
		print,'mean=',mean
;		yr=yr-mean
                ymin=min([yr,yi])
                ymax=max([yr,yi])
                set_xy,xmin,xmax,ymin,ymax
                plot,xpow,yr
                !p.color=truecolor(ncolor*3/4)
                oplot,xpow,yi

; panel 4: power spectrum: phi=exp[-i(omega*t+m*theta-n*zeta)]
        !noeras=1
       set_viewport,0.59,0.99,0.05,0.45
                !p.color=truecolor(ncolor+1)
		power=complex(yr,yi)
                power=fft(power,-1)
                ypow=abs(power)
                yp=fltarr(2*np)
                xp=fltarr(2*np)
                for i=0,np-2 do begin
                        yp(i)=ypow(i+ntime-nbegin-np+1)
                        xp(i)=(i-np+1)*6.283185/((ntime-nbegin)*delt)
                end
                for i=0, np do begin
                        yp(np-1+i)=ypow(i)
                        xp(np-1+i)=i*6.283185/((ntime-nbegin)*delt)
                end

	m=0
;	print,'m=?'
;	read,m
	ktheta=m/120.0
	shift=ktheta*0.005
	print,shift
	xp=xp+shift

                !mtitle="frequency spectrum"
                xmax=max(xp)
                xmin=min(xp)
                ymax=max(yp)
                ymin=min(yp)
;		xmax=0.008
;		xmin=-0.008
                set_xy,xmin,xmax,ymin,ymax
                plot,xp,yp,/xstyle

        if !D.name eq 'PS' then begin
                device,/close
                set_plot,'X'
		truecolor=true_color
		!p.thick=1
        endif

return
end
;****************************************************************
pro hplt

common startup,number_plot,fpath,color_max,color_value,direct_color
common idlh,plotidh,ncolor,np,nm,glength,truecolor,pseudo_color,true_color
common data,f,hname,xtime,nbegin,naverage,ntime,delt

plotidh=number_plot
ncolor=color_max
true_color=color_value
truecolor=true_color
pseudo_color=direct_color

; default window
        set_plot,'X'
	!p.thick=1
        !linetype=0
        !p.background=truecolor(ncolor/2)
        !p.color=truecolor(ncolor+1)
        window, plotidh, TITLE='history', xsize=600,ysize=600

	fname = dialog_pickfile(/read, path=fpath, filter='*.out')
	if fname eq '' then begin
	wdelete
	return
	endif
	openr, plotidh, fname

;openr, plotidh, 'history.out'

	nmode=5
        nradial=5
	nother=18
	nrun=0
	ntime=1
	nbegin=0
	naverage=0
	delt=1.0
; nplot=# of plots, nradial=# of radial bin, ntime=# of time step
; nlast=int. points of flux surface

        readf,plotidh, nrun,nother,nradial,nmode,ntime,delt

	print,'total time step=',ntime
        print,'read # of last time step'
	new=''
        read,new
	if new ne '' then ntime=0+new
;	nbegin=1
;	print,"read # of first time step"
;	new=''
;        read,new
;	if new ne '' then nbegin=0+new

	f=fltarr(4*nmode+nradial+nother,ntime)
	readf, plotidh, f
	close, plotidh

	rmajor=1.0/0.00142
	xtime=indgen(ntime)
	np=(ntime-nbegin)/50
	if np lt 2 then np=2

	nm=(ntime-nbegin)/4
	glength=1.0

;        glength=527.0
;        glength=404.0

hname=strarr(45)
hname=["Exit","???","orbit poloidal", "orbit flux", "weight",$
	"V_para","energy","angular momentum", "ddeni","ddene",$
	"radial mode", "field energy", "entropyi","entropye","dflowi",$
	"dlfowe","pfluxi","pfluxe","efluxi","efluxe",$
	"local fluxes","shear rate1","shear rate2","shear rate3",$
	"shear rate4","shear rate5","shear rate6","shear rate7",$
	"shear rate8","(m,n) mode1 ","(m,n) mode2","(m,n) mode3","(m,n) mode4",$
	"(m,n) mode5","(m,n) mode6","(m,n) mode7","(m,n) mode8","peak chi_i",$
	"nbegin,naverage","frequency","X window2","X_size","PS file"]

xmenu,hname,BASE=hbase,SPACE=5,TITLE='history',column=2,xpad=5,ypad=5
widget_control,hbase,/realize
xmanager,"history",hbase,/NO_BLOCK

end
;**************************************************************************
