pro snapshot_event,event

common idls,plotids,ncolor,layer,cvalue,mode0,cutoff,nkr,nkz,truecolor,pseudo_color,true_color
common sdata,step,q,mbin_u,mbin_psi,mpsi,mtheta,mzetamax,nflux,data_u,$
	data_r,data_p,data_f,m_toroidal,m_poloidal,mode_eigen,data_eigen,position

widget_control,event.id,get_uvalue=choice

case choice of
 
0: begin            ; Quit
	set_plot,'X'
        wdelete        ; Closes plot windowse
        widget_control, event.top,/destroy
end

1: begin		;f(u_para)	
	title=['f(u_para)','delta f(u_para)']
	x=data_u(*,0)
	y1=data_u(*,(position-1)*6+3)
	y2=data_u(*,(position-1)*6+4)
	plot2frame,x,y1,x,y2,title
end	

2: begin		;f(pitch)
	title=['f(pitch)','delta f(pitch)']
	x=data_u(*,1)
	y1=data_u(*,(position-1)*6+5)
	y2=data_u(*,(position-1)*6+6)
	plot2frame,x,y1,x,y2,title
end	

3: begin		;f(energy)	
	title=['f(energy)','delta f(energy)']
	x=data_u(*,2)
	y1=data_u(*,(position-1)*6+7)
	y2=data_u(*,(position-1)*6+8)
	plot2frame,x,y1,x,y2,title
	print,position
end	

4: begin		;density	
	title=['zonali','zonale']
	x=data_r(*,0)
	y1=data_r(*,1)
	y2=data_r(*,2)
	plot2frame,x,y1,x,y2,title	
end	

5: begin		;sheared flow	
	title=['zonal flow','marker']
	x=data_r(*,0)
	y1=data_r(*,3)
	y2=data_r(*,4)
	
	plot2frame,x,y1,x,y2,title
end	

6: begin		;parallel	
	title=['u_para(r)', 'delta u_para(r)']
	x=data_r(*,0)
	y1=data_r(*,5)
	y2=data_r(*,6)
	plot2frame,x,y1,x,y2,title
end	

7: begin		;temperature	
	title=['T(r)', 'delta T(r)']
	x=data_r(*,0)
	y1=data_r(*,7)
	y2=data_r(*,8)
	plot2frame,x,y1,x,y2,title
end	

8: begin		;marker
	!mtitle=''	
	x=data_p(*,*,0)
	y=data_p(*,*,1)
	f=data_p(*,*,2)
	poloidal,x,y,f,mpsi,mtheta
end	

9: begin		;density	
	!mtitle=''	
	x=data_p(*,*,0)
	y=data_p(*,*,1)
	f=data_p(*,*,2)
	poloidal,x,y,f,mpsi,mtheta
end	

10: begin		;potential	
;	!mtitle='phi(r,theta)'	
	x=data_p(*,*,0)
	y=data_p(*,*,1)
	f=data_p(*,*,2)
;	for j=0,mtheta do begin
;		print,j,f(mpsi-1,j)
;	endfor
		
	poloidal,x,y,f,mpsi,mtheta
end
	
11: begin		;phi(r)	
	title=['phi(r)', 'sqrt(phi^2)(r)']	
	f=data_p(*,*,2)
	x=data_r(*,0)
	y1=fltarr(mpsi)
	y2=fltarr(mpsi)
        for i=0,mpsi-1 do begin
                y1(i)=0.0
		y2(i)=0.0
                for j=0,mtheta do begin
                        y1(i)=y1(i)+f(i,j)
                        y2(i)=y2(i)+f(i,j)*f(i,j)
                endfor
                y1(i)=y1(i)/mtheta
                y2(i)=sqrt(y2(i)/mtheta)
        endfor
	plot2frame,x,y1,x,y2,title
end	

12: begin		;phi(theta)	
	title=['phi(theta)', 'sqrt(phi^2)(theta)']
	f=data_p(*,*,2)
        x=(indgen(mtheta)/float(mtheta-1)-0.5)*2.0*3.14
        y1=fltarr(mtheta)
        y2=fltarr(mtheta)
        for j=0,mtheta-1 do begin
                y1(j)=f(nflux,j)
		y2(j)=0.0
		for i=0,mpsi-1 do begin	
			y2(j)=y2(j)+f(i,j)*f(i,j)
		endfor
		y2(j)=sqrt(y2(j)/float(mpsi))
        endfor
	plot2frame,x,y1,x,y2,title
end	

13: begin		;k_r	
	title=['k_r', 'k_r']
	f=data_p(*,*,2)
        x1=indgen(nkr)
        y1=fltarr(nkr)*0.0
	x2=data_r(*,0)
	y2=x2
	for i=1,mpsi-1 do begin
		y2(i)=x2(i)-x2(i-1)
	endfor
	y2(0)=y2(1)

        yy=fltarr(mpsi)
        for j=mtheta/4,3*mtheta/4 do begin
		yy=f(*,j)
                yy=fft(yy,1)
               	y1(0)=y1(0)+abs(yy(0))
               	for i=1,nkr-1 do begin
                       	y1(i)=y1(i)+sqrt(yy(i)^2+yy(mpsi-i)^2)
               	endfor
        endfor

	sum=0.5*(y1(0)+y1(nkr-1))
	for i=1,nkr-1 do begin
		sum=sum+y1(i)
	endfor
	y1=y1*(nkr-1)/sum

;	krave=0.0
;	insten=0.0
;	for i=0,16 do begin
;		krave=krave+i*i*y1(i)*y1(i)
;		insten=insten+y1(i)*y1(i)
;	endfor
;	print,sqrt(krave/insten)
	
	krdat='kr.dat'
	kr=''
;	print, 'file name=?'
;	read, kr
	if kr ne '' then krdat=kr
;	openw,111,krdat
;	printf, 111, y1
;	close,111

	plot2frame,x1,y1,x2,y2,title
end	

14: begin		;phi(zeta,theta)	
;	!mtitle='phi(zeta,theta)'
	x=2.0*3.1415926*((indgen(mzetamax)+1.0)/mzetamax-0.5)
	y=2.0*3.1415926*((indgen(mtheta)+1.0)/mtheta-0.5)
	f=data_f(0:mzetamax-1,0:mtheta-1,0)

	ff=f
	for i=0,mzetamax-1 do begin
		ii=i-mzetamax/2
		if ii lt 0 then ii=ii+mzetamax
;		ff(ii,*)=f(i,*)
	endfor
	f=ff

		for j=0,mtheta-1 do begin
			jj=j-mtheta/2
			if jj lt 0 then jj=jj+mtheta
;			ff(*,jj)=f(*,j)
		endfor

	flux,x,y,ff,mzetamax,mtheta
end	

15: begin		;k_theta,k_zeta	
	title=['k_theta','k_zeta']
	f=data_f(0:mzetamax-1,0:mtheta-1,0)

        x1=indgen(nkz)
        y1=fltarr(nkz)*0.0
; k_theta
        yy=fltarr(mtheta)
        for i=0,mzetamax-1 do begin
        	yy=f(i,*)
                yy=fft(yy,1)
                y1(0)=y1(0)+(abs(yy(0)))^2
              	for j=1,nkz-1 do begin
                       	y1(j)=y1(j)+(abs(yy(j)))^2+(abs(yy(mtheta-j)))^2
                endfor
        endfor
	y1=sqrt(y1/mzetamax)
	y1=y1/mtheta
;	print,y1

; k_para
        x2=indgen(8)
        y2=fltarr(8)*0.0
	mconn=fix(float(mzetamax)*q+0.5)
        yy=fltarr(mconn)
        for i=0,mtheta-1 do begin
		for j=0,mconn-1 do begin
			it=i
			jt=j
			if jt gt mzetamax-1 then begin
				jt=jt-mzetamax
				it=it+fix(float(mtheta)/q+0.5)
				if it gt mtheta-1 then it=it-mtheta
			endif
	        	yy(j)=f(jt,it)
		endfor
                yy=fft(yy,1)
                y2(0)=y2(0)+(abs(yy(0)))^2
              	for j=1,7 do begin
                       	y2(j)=y2(j)+(abs(yy(j)))^2+(abs(yy(mconn-j)))^2
                endfor
        endfor
	y2=sqrt(y2/mtheta)
	y2=y2/mconn	
;	print,y2

	plot2frame,x1,y1,x2,y2,title
;	plot2frame,x1(300:320),y1(300:320),x2,y2,title
end	

16: begin		;k_theta,k_para	
; for non-field-line-following
;	title=['k_theta', 'k_para']	
	f=data_f(0:mzetamax-1,0:mtheta-1,0)

        x1=indgen(nkz)
        y1=fltarr(nkz)*0.0
; k_theta
        yy=fltarr(mtheta)*0.0
        for i=0,mzetamax-1 do begin
        	yy=f(i,*)
	        yy=fft(yy,1)
	        y1(0)=y1(0)+(abs(yy(0)))^2
	        for j=1,nkz-1 do begin
	        	y1(j)=y1(j)+(abs(yy(j)))^2+(abs(yy(mtheta-j)))^2
	        endfor
	endfor
	y1=sqrt(y1/mzetamax)
	y1=y1/mzetamax

; k_para
        x2=indgen(nkz/4)
        y2=fltarr(nkz/4)*0.0
; use linear interpretation
	mconn=fix(float(mzetamax)*q+0.5)
        yy=fltarr(mconn)*0.0
	for jt=0,mtheta-1 do begin
		for it=0,mconn-1 do begin
			ratio=float(mtheta)/(q*float(mzetamax))
			j=jt+fix(float(it)*ratio)
			delt=float(jt)+float(it)*ratio-float(j)
			jp=j+1
			if j gt mtheta-1 then j=j-mtheta
			if jp gt mtheta-1 then jp=jp-mtheta
			i=it
			if i gt mzetamax-1 then i=i-mzetamax

			yy(it)=delt*f(i,jp)+(1.0-delt)*f(i,j)
		endfor
		yy=fft(yy,1)
		y2(0)=y2(0)+(abs(yy(0)))^2
	        for j=1,nkz/4-1 do begin
        		y2(j)=y2(j)+(abs(yy(j)))^2+(abs(yy(mconn-j)))^2
	        endfor
        endfor
	y2=sqrt(y2/mtheta)
	y2=y2/mconn

	plot2frame,x1,y1,x2,y2,title

;        x=1.0*indgen(mzetamax)
;        y1=x
;	for i=0,mzetamax-1 do begin
;		y1(i)=f(i,0)*f(i,0)
;		for j=1,mtheta-1 do begin
;			y1(i)=y1(i)+f(i,j)*f(i,j)		
;		endfor
;	endfor
;	y1=sqrt(y1/mtheta) 
;	y2=y1
;	plot2frame,x,y1,x,y2,title

end	

17: begin	;1st mode
	m=0
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

18: begin	;2nd mode
	m=1
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

19: begin	;3rd mode
	m=2
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

20: begin	;4th mode
	m=3
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

21: begin	;5th mode
	m=4
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

22: begin	;6th mode
	m=5
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

23: begin	;7th mode
	m=6
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

24: begin	;8th mode
	m=7
	x=data_r(*,0)
	plot_eigen,m,mode_eigen(m),mpsi,m_poloidal,x,data_eigen
end	

25: begin	;plot to PS
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

26: begin	;subtract (0,1) mode
	print, 'old mode0=',mode0
	print, 'new mode0=?'
	read, mode0	
end	

27: begin	;cutoff in contour plot
	print, 'old cutoff=',cutoff
	print, 'new cutoff=?'
	read, cutoff	
end	

28: begin	;change X-window size
	size_x=1
	size_y=1
	print, 'size_x,size_y=?'
	read,size_x,size_y
        window, plotids, title='snapshot'+string(plotids), xsize=size_x,ysize=size_y
end	

29: begin	;change number of layer for contour plot
;	print, 'old layer=',layer
;	print, 'new layer=?'
;	read, layer
	read, position
;	cvalue=intarr(layer)
;	cvalue=indgen(layer)*ncolor/(layer-1)
;	print,'color table=',cvalue
end	

30: begin	;change number of kr spectral
	print, 'old nkr=',nkr
	print, 'new nkr=?'
	read, nkr
end	

31: begin	;change number of kz spectral
	print, 'old nkz=',nkz
	print, 'new nkz=?'
	read, nkz
end	

endcase
end
;****************************************************************
pro plot_eigen,m,mode,mpsi,m_poloidal,x,data_eigen

common idls,plotids,ncolor,layer,cvalue,mode0,cutoff,nkr,nkz,truecolor,pseudo_color,true_color

	if !D.name eq 'X' then wset,plotids
	!noeras=0
	!p.background=truecolor(ncolor/2)
	!p.color=truecolor(ncolor+1)

	xmax=max(x)
        xmin=min(x)
;	xmin=20
;	xmax=80
	ymax=max(data_eigen(*,m,*))
        ymin=min(data_eigen(*,m,*))
        set_xy,xmin,xmax,ymin,ymax	
	set_viewport,0.1,0.9,0.1,0.9
	!mtitle='n='+string(mode)
	
	ltype=[0,5,0,2,0,2,0,5,0]
	lcolor=[ncolor*5/16,ncolor/16,ncolor/8,ncolor/4,ncolor+1,ncolor*3/4-1,ncolor*7/8-1,ncolor*15/16-1,ncolor*11/16-1]

	y1=fltarr(mpsi)
	for i=0,m_poloidal-1 do begin
		y1=data_eigen(i,m,*)

		!linetype=ltype(i)
		!p.color=truecolor(lcolor(i))
;		print, lcolor(i),truecolor(lcolor(i))
		if i eq 0 then plot,x,y1,/xstyle
		if i gt 0 then oplot,x,y1
	endfor

	if !D.name eq 'PS' then begin
		device,/close
		set_plot,'X'
		truecolor=true_color
	        !p.thick=1
	endif

return
end
;****************************************************************
pro plot2frame,x1,y1,x2,y2,title
common idls,plotids,ncolor,layer,cvalue,mode0,cutoff,nkr,nkz,truecolor,pseudo_color,true_color

	if !D.name eq 'X' then wset,plotids
	!noeras=0
	!p.background=truecolor(ncolor/2)
	!p.color=truecolor(ncolor+1)
	!linetype=0
	
	set_viewport,0.14,0.54,0.1,0.9
;	set_viewport,0.1,0.9,0.1,0.9
	xmax=max(x1)
        xmin=min(x1)
	ymax=max(y1)
        ymin=min(y1)

;	ymax=0.01
;	ymin=0.000001
;	xmin=0.0
;	xmax=160.0
        set_xy,xmin,xmax,ymin,ymax	
	!mtitle=title(0)
;	print,y1
	plot,x1,y1,/xstyle
;	!linetype=1
;	y1=y1(0)*exp(-x1/6.0)
;	oplot,x1,y1
;return

	!linetype=0
	!noeras=1
	set_viewport,0.59,0.99,0.1,0.9
;	set_viewport,0.09,0.99,0.1,0.9
;	y2=alog(y2)
;	x2=alog(x2)
;	x2(0)=1
	xmax=max(x2)
        xmin=min(x2)
	ymax=max(y2)
        ymin=min(y2)

;	xmin=0.0
;	xmax=160.0
;	ymax=0.01
;	ymin=0.0
        set_xy,xmin,xmax,ymin,ymax	
	!mtitle=title(1)
;	print,y2
	plot,x2,y2,/xstyle
;	plot,/ylog,/xlog,x2,y2

	if !D.name eq 'PS' then begin
		device,/close
		set_plot,'X'
		truecolor=true_color
	        !p.thick=1
	endif

;        openw,1,'data.dat'
;        printf,1,x1,y1
        close,1

return
end
;****************************************************************
pro poloidal,x,y,f,mpsi,mtheta

common idls,plotids,ncolor,layer,cvalue,mode0,cutoff,nkr,nkz,truecolor,pseudo_color,true_color

	if max(f) eq min(f) then print,max(f)
	if max(f) eq min(f) then return

	if !D.name eq 'X' then wset,plotids
	!noeras=0
	!p.background=truecolor(ncolor/2)
	!p.color=truecolor(ncolor+1)
; 8-bits color
;        !p.background=ncolor
;        !p.color=ncolor+1

; cover axis
;	x(0,0:mtheta)=1.0
;	y(0,0:mtheta)=0.0
	f(0,0:mtheta)=0.0
	f(mpsi-1,0:mtheta)=0.0

; substract (0,1) mode
	pi=3.1415927
	mode01=0.0*indgen(mpsi)
	mode10=0.0*indgen(mpsi)
	mode00=0.0*indgen(mpsi)
;        for j=0,mtheta-1 do begin
;            tdum=2.0*pi*(j+0.5)/mtheta
;            tsin=sin(tdum)
;            tcos=cos(tdum)
;            for i=0,mpsi-1 do begin
;                mode00(i)=mode00(i)+f(i,j)
;                mode01(i)=mode01(i)+tcos*f(i,j)
;;                mode10(i)=mode10(i)+tsin*f(i,j)
;            endfor
;        endfor
;        for i=0,mpsi-1 do begin
;            mode00(i)=mode00(i)/mtheta
;            mode01(i)=mode01(i)*2.0/mtheta
;            mode10(i)=mode10(i)*2.0/mtheta
;            mode00(i)=mode00(i)/mtheta
;        endfor
;        xmax=0
;        xmin=mpsi-1
;        ymax=max(mode00)
;        ymin=min(mode00)
;        set_xy,xmin,xmax,ymin,ymax
;        plot,mode00
;        read,dump

;        xmax=0
;        xmin=mpsi-1
;        ymax=max(mode00)
;        ymin=min(mode00)
;        set_xy,xmin,xmax,ymin,ymax
;        plot,mode01
;        read,dump

;        xmax=0
;        xmin=mpsi-1
;        ymax=max(mode00)
;        ymin=min(mode00)
;        set_xy,xmin,xmax,ymin,ymax
;        plot,mode10
;        read,dump

;	for i=0,mpsi-1 do begin
;		for j=0,mtheta do begin
;                    tdum=2.0*pi*(j+0.5)/mtheta
;                    f(i,j)=f(i,j)-mode00(i)
;                    f(i,j)=f(i,j)-mode01(i)*cos(tdum)*(1.0-mode0)
;                    f(i,j)=f(i,j)-mode10(i)*sin(tdum)*(1.0-mode0)
;		endfor
;	endfor

;	fp=fltarr(mpsi,mtheta+1)
;	fm=fltarr(mpsi,mtheta+1)

;for i=1,2 do begin
; radial smoothing
;	fp(0:mpsi-2,0:mtheta)=f(1:mpsi-1,0:mtheta)
;	fp(mpsi-1,0:mtheta)=f(mpsi-1,0:mtheta)
;	fm(1:mpsi-1,0:mtheta)=f(0:mpsi-2,0:mtheta)
;	fm(0,0:mtheta)=f(0,0:mtheta)
;	f=0.5*f+0.25*(fp+fm)

; poloidal smoothing
;	fp(0:mpsi-1,0:mtheta-1)=f(0:mpsi-1,1:mtheta)
;	fp(0:mpsi-1,mtheta)=f(0:mpsi-1,0)
;	fm(0:mpsi-1,1:mtheta)=f(0:mpsi-1,0:mtheta-1)
;	fm(0:mpsi-1,0)=f(0:mpsi-1,mtheta)
;	f=0.5*f+0.25*(fp+fm)
;endfor	

; f value	
;	for i=0,mpsi-1 do begin
;		for j=0,mtheta do begin
;			if f(i,j) gt 0.8 then f(i,j)=0.0
;		endfor
;	endfor

; select poloidal spectrum
        jmax=mtheta/4
        for i=0,mpsi-1 do begin
            pd=f(i,0:mtheta-1)
            pf=float(mtheta)*fft(pd)
            pf(jmax:mtheta-1-jmax)=0.0
            pd=fft(pf,-1)
            for j=0,mtheta-1 do begin
                f(i,j)=pd(mtheta-1-j)
            endfor
        endfor
        f(0:mpsi-1,mtheta)=f(0:mpsi-1,0)

        imax=mpsi/4
        for j=0,mtheta do begin
            rd=f(0:mpsi-1,j)
            rf=float(mpsi)*fft(rd)
            rf(imax:mpsi-1-imax)=0.0
            rd=fft(rf,-1)
            for i=0,mpsi-1 do begin
                f(i,j)=rd(mpsi-1-i)
            endfor
        endfor

	zmax=max(f)
	zmin=min(f)
        print,zmin,zmax
	flevel=zmin+(zmax-zmin)*indgen(layer)/(layer-1)

; ploting range
        length=2793.3*2.0
;        length=1.0
        x=x*length
        y=y*length
        xmax=max(x)
        xmin=min(x)
        ymax=1500.0
        ymin=-1500.0
;        del=(ymax-ymin-xmax+xmin)/2.0
;	del=0.0
;        if del gt 0.0 then xmax=xmax+0.8*del
;        if del gt 0.0 then xmin=xmin-0.8*del
;        if del lt 0.0 then ymax=ymax-0.8*del
;        if del lt 0.0 then ymin=ymin+0.8*del
        set_xy,xmin,xmax,ymin,ymax

; surface
;        set_viewport,0.08,0.48,0.1,0.9
;        surface,f,x,y,az=30,charsize=2.0
;        !noeras=1
; contour
;        set_viewport,0.16,0.96,0.1,0.9
        set_viewport,0.075,0.985,0.055,0.965

;	if !D.name eq 'X' then 
;	polyfill, [0,1,1,0,0],[0,0,1,1,0],/normal, color=ncolor
;        !noeras=1
	!x.style=1
	!y.style=1
        contour,f,x,y,nlevels=layer,c_colors=cvalue,max_value=zmax,min_value=zmin,$
		levels=flevel,/fill,charsize=1
 
;openr,1,'orbit.dat'
;nt=1
;readf,1,nt
;print,nt
;xt=fltarr(nt)
;yt=fltarr(nt)
;readf,1,xt,yt
;xt=(xt+1.0)*length
;yt=yt*length
;oplot,xt,yt
;close,1

;        !noeras=1
;	oplot,x(mpsi-1,0:mtheta),y(mpsi-1,0:mtheta)

	if !D.name eq 'PS' then begin
		device,/close
		set_plot,'X'
		truecolor=true_color
	        !p.thick=1
	endif

return
end
;****************************************************************
pro flux,x,y,f,mzetamax,mtheta

common idls,plotids,ncolor,layer,cvalue,mode0,cutoff,nkr,nkz,truecolor,psedocolor,true_color

	if !D.name eq 'X' then wset,plotids
	!noeras=0
	!p.background=truecolor(ncolor/2)
	!p.color=truecolor(ncolor+1)

;        xmax=max(x)
;        xmin=min(x)
;        ymax=max(y)
;        ymin=min(y)
;        del=(ymax-ymin-xmax+xmin)/2.0
;        if del gt 0.0 then xmax=xmax+0.8*del
;        if del gt 0.0 then xmin=xmin-0.8*del
;        if del lt 0.0 then ymax=ymax-0.8*del
;        if del lt 0.0 then ymin=ymin+0.8*del
	xmin=-3.5
	xmax=3.5
	ymin=-3.5
	ymax=3.5
        set_xy,xmin,xmax,ymin,ymax

; substract (0,1) mode
;	pi=3.1415927
;	mode01=0.0
;	mode10=0.0
;	for j=0,mtheta-1 do begin
;		tdum=cos(y(j))
;		tdum1=sin(y(j))
;		for i=0,mzetamax-1 do begin
;			mode01=mode01+tdum*f(i,j)
;			mode10=mode10+tdum1*f(i,j)
;		endfor
;	endfor
;	mode01=mode01*2.0/mzetamax/mtheta
;	for j=0,mtheta-1 do begin
;		tdum=mode01*cos(y(j))
;		tdum1=mode10*sin(y(j))
;		tdum1=0.0
;		for i=0,mzetamax-1 do begin
;			f(i,j)=f(i,j)-(tdum+tdum1)*(1.0-mode0)
;		endfor
;	endfor

;	fp=fltarr(mzetamax,mtheta)
;	fm=fltarr(mzetamax,mtheta)
; toroidal smoothing
;	fp(0:mzetamax-2,0:mtheta-1)=f(1:mzetamax-1,0:mtheta-1)
;	fp(mzetamax-1,0:mtheta-1)=f(0,0:mtheta-1)
;	fm(1:mzetamax-1,0:mtheta-1)=f(0:mzetamax-2,0:mtheta-1)
;	fm(0,0:mtheta-1)=f(mzetamax-1,0:mtheta-1)
;	f=0.5*f+0.25*(fp+fm)

; poloidal smoothing
;	fp(0:mzetamax-1,0:mtheta-2)=f(0:mzetamax-1,1:mtheta-1)
;	fp(0:mzetamax-1,mtheta-1)=f(0:mzetamax-1,0)
;	fm(0:mzetamax-1,1:mtheta-1)=f(0:mzetamax-1,0:mtheta-2)
;	fm(0:mzetamax-1,0)=f(0:mzetamax-1,mtheta-1)
;	f=0.5*f+0.25*(fp+fm)

; f range	
	zmax=max(f)
	zmin=min(f)
	zmax=(0.5-cutoff)*(zmax-zmin)
	zmin=-zmax

; surface
;        set_viewport,0.2,0.5,0.1,0.9
;        surface,f,x,y,az=30,charsize=2.0
;        !noeras=1
; contour

        set_viewport,0.1,0.9,0.1,0.9
        contour,f,x,y,nlevels=layer,c_colors=cvalue,max_value=zmax,min_value=zmin,/fill,/xstyle,/ystyle

	if !D.name eq 'PS' then begin
		device,/close
		set_plot,'X'
		truecolor=true_color
	        !p.thick=1
	endif

return
end
;****************************************************************
pro splt

common startup, number_plot,fpath,color_max,color_value,direct_color
common idls,plotids,ncolor,layer,cvalue,mode0,cutoff,nkr,nkz,truecolor,pseudo_color,true_color
common sdata,step,q,mbin_u,mbin_psi,mpsi,mtheta,mzetamax,nflux,data_u,$
	data_r,data_p,data_f,m_toroidal,m_poloidal,mode_eigen,data_eigen,position

true_color=color_value
truecolor=true_color
pseudo_color=direct_color
plotids=number_plot
ncolor=color_max
layer=color_max
cvalue=intarr(layer)
cvalue=truecolor(indgen(layer)*ncolor/layer)
;8-bits color
;cvalue=indgen(layer)*ncolor/(layer-1)
mode0=0.0
cutoff=-1.0
position=1

; default window
	set_plot,'X'
	!linetype=0
	!p.thick=1
	!p.background=truecolor(ncolor/2)
	!p.color=truecolor(ncolor+1)
        window,  plotids, title='snapshot'+string(plotids), xsize=1000,ysize=1000

fname = dialog_pickfile(/read, path=fpath, filter='*.out')
if fname eq '' then begin
	wdelete
	return
endif
openr, plotids, fname

; # of snapshots and cell dimension
	m=1
	mbin_u=1
	mbin_psi=1
	mpsi=1
	mtheta=1
	mzetamax=1
	nflux=1
	readf, plotids, step,q,mbin_u,mbin_psi,mpsi,mtheta,mzetamax,nflux

	nkr=mpsi/6
	nkz=mtheta/4

; velocity space distribution
	readf, plotids, m
	data_u=fltarr(mbin_u,mbin_psi*(m-3)+3)
	readf, plotids, data_u

; radial profile
	readf, plotids, m
	data_r=fltarr(mpsi,m)
	readf, plotids, data_r

; poloidal cross section profile
	readf, plotids, m
	print,mpsi,mtheta,m
	data_p=fltarr(mpsi,mtheta+1,m)
	readf, plotids, data_p

; flux surface contour plot
	readf, plotids, m
	data_f=fltarr(mzetamax,mtheta,m)
	readf, plotids, data_f

; eigenmode structure
	m_poloidal=1
	mmode=intarr(mpsi)
;	readf, plotids, m_poloidal,mmode
;	print,'toroidal mode #=',mmode

	data_eigen=fltarr(m_poloidal,2,mpsi)
;	readf, plotids, data_eigen
	close, plotids

; make a menu
sname=strarr(32)
sname=["Exit","f(v_para","f(pitch)","f(energy)","density(r)",$
	"shear flow","u_para(r)","Temperature(r)","n(r,theta)","delta_n(r,theta)",$
	"phi(r,theta)","phi(r)","phi(theta)","k_r","phi(theta,zeta)",$
	"k_theta,k_zeta","k_theta,k_para","mode1","mode2","mode3",$
	"mode4","mode5","mode6","mode7","mode8",$
	"PS file","mode01","cutoff","X_size","# of layer","nkr","nkz"]

xmenu,sname,BASE=sbase,SPACE=5,title='snapshot'+string(plotids),column=2,xpad=5,ypad=5
widget_control,sbase,/realize
xmanager,"snapshot",sbase,/NO_BLOCK

end
;********************************************************************************




