pro equilibrium_event,event

common idle,plotide,ncolor,layer,cvalue,truecolor,pseudo_color,true_color,nind0,nind1
common edata,ename,nrplot,nplot,lsp1,lsp2,lst,isp,pdata,tdata,spdata


widget_control,event.id,get_uvalue=choice

case choice of

0: begin		;minor radius
    title='inverse aspec-ratio from profile data'
    x=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,2)
    y2=y1	
    panel1,x,y1,x,y2,title
end	

1: begin		;major radius
    title='major radius from profile data'
    x=pdata(nind0:nind1,0)

    y1=pdata(nind0:nind1,3)
    y2=y1
    panel1,x,y1,x,y2,title
end	

2: begin		;Te
    title1='-dlnTe/dpsi'
    title2='Te(psi)'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23)
    y1=pdata(nind0:nind1,5)
    x2=pdata(nind0:nind1,0)
    y2=pdata(nind0:nind1,4)
    panel2,x1,y1,x2,y2,title1,title2
end	

3: begin		;ne
    title1='-dlnne/dpsi'
    title2='ne(psi)'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23)
    y1=pdata(nind0:nind1,7)
    x2=pdata(nind0:nind1,0)
    y2=pdata(nind0:nind1,6)
    panel2,x1,y1,x2,y2,title1,title2
end	

4: begin		;Ti
    title1='-dlnTi/dpsi'
    title2='Ti(psi)'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23)
    y1=pdata(nind0:nind1,9)
    x2=pdata(nind0:nind1,0)
    y2=pdata(nind0:nind1,8)
    panel2,x1,y1,x2,y2,title1,title2
end	

5: begin		;ni
    title1='-dlnni/dpsi'
    title2='ni(psi)'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23);
    y1=pdata(nind0:nind1,11);
    x2=pdata(nind0:nind1,0);
    y2=pdata(nind0:nind1,10);
    panel2,x1,y1,x2,y2,title1,title2
    print,'avg dlnni/dpsi=',mean(y1)
end	

6: begin		;Tf
    title1='-dlnTf/dpsi'
    title2='Tf(psi)'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23)
    y1=pdata(nind0:nind1,13)
    x2=pdata(nind0:nind1,0)
    y2=pdata(nind0:nind1,12)
    panel2,x1,y1,x2,y2,title1,title2
end	

7: begin		;nf
    title1='-dlnnf/dpsi'
    title2='nf(psi)'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23)
    y1=pdata(nind0:nind1,15)
    x2=pdata(nind0:nind1,0)
    y2=pdata(nind0:nind1,14)
    panel2,x1,y1,x2,y2,title1,title2
end

8: begin		;Zeff
    title=ename(choice)
    x=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,16)
    y2=y1
    panel1,x,y1,x,y2,title
end	

9: begin		;rotation
    title=ename(choice)
    x=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,17)
    y2=y1
    panel1,x,y1,x,y2,title
end	

10: begin		;E_r
    title=ename(choice)
    x=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,18)
    y2=y1
    panel1,x,y1,x,y2,title
end	

11: begin		;q(psi),q(r)
    title1='q(psi)'
    title2='q(r)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,19)
    x2=pdata(nind0:nind1,23)
    y2=y1
    panel2,x1,y1,x2,y2,title1,title2
end	

12: begin		;shear(psi),s(r)
    title1='psi(r)'
    title2='dq/dpsi'
    x1=pdata(nind0:nind1,23)/pdata(lsp1-1,23)    
    x2=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,0)
    y2=pdata(nind0:nind1,20)
    panel2,x1,y1,x2,y2,title1,title2
end

13: begin		;current g(psi) & g(r)
    title1='g(psi)'
    title2='g(r)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,21)
    x2=pdata(nind0:nind1,23)
    y2=y1
    panel2,x1,y1,x2,y2,title1,title2
end	

14: begin		;pressure p(psi) & p(r)
    title1='P(psi)'
    title2='P(r)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,22)
    x2=pdata(nind0:nind1,23)
    y2=y1
    panel2,x1,y1,x2,y2,title1,title2
end	

15: begin		;toroidal flux (psi) & (r)
    title1='toroidal flux(psi)'
    title2='toroidal flux(r)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,24)
    x2=pdata(nind0:nind1,23)
    y2=y1
    panel2,x1,y1,x2,y2,title1,title2
end	

16: begin		;r(psi) & r(torpsi)
    title1='r(psi)'
    title2='r(toroidal flux)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,23)
    x2=pdata(nind0:nind1,24)
    y2=y1
    panel2,x1,y1,x2,y2,title1,title2
end	

17: begin		;radial grid
    title1='radial grid rg(psi)'
    title2='radial grid rg(r)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,25)
    x2=pdata(nind0:nind1,23)
    y2=y1
    panel2,x1,y1,x2,y2,title1,title2
end	

18: begin		;toroidal flux and its inverse
    title1='toroidal flux tor(psi) & psi(tor)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,24)
    x2=pdata(0,24)+(pdata(lsp1-1,24)-pdata(0,24))*indgen(lsp1)/(lsp1-1.0)
    y2=pdata(nind0:nind1,26)
    panel1,x1,y1,y2,x2,title1
end	

19: begin		;radial grid and its inverse
    title1='radial grip rg(psi) & psi(rg)'
    x1=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,25)
    x2=pdata(0,25)+(pdata(lsp1-1,25)-pdata(0,25))*indgen(lsp1)/(lsp1-1.0)
    y2=pdata(nind0:nind1,27)
    panel1,x1,y1,y2,x2,title1
end

20: begin		;q(psi)
    title='q(psi) from spdata & dtorpsi/dpsi'
    x=pdata(nind0:nind1,0)
    y1=pdata(nind0:nind1,19)
    y2=fltarr(lsp1)
    for i=1,lsp1-2 do begin
        y2(i)=(pdata(i+1,24)-pdata(i-1,24))/(pdata(i+1,0)-pdata(i-1,0))
    endfor
    y2(0)=(pdata(1,24)-pdata(0,24))/(pdata(1,0)-pdata(0,0))
    y2(lsp1-1)=(pdata(lsp1-1,24)-pdata(lsp1-2,24))/(pdata(lsp1-1,0)-pdata(lsp1-2,0))
    panel1,x,y1,x,y2,title
end

21: begin		;psi(rg) & r(rg)
    title1='radial grid psi(rg)'
    title2='radial grid dpsi/dr'
    x1=pdata(0,25)+(pdata(lsp1-1,25)-pdata(0,25))*indgen(lsp1)/(lsp1-1.0)
    y1=pdata(nind0:nind1,27)
    x2=x1
    y2=DERIV(x1,y1)
    ddy=mean(DERIV(x1,y2))
    panel2,x1,y1,x2,y2,title1,title2
    print,"For c*r=dpsi/dr,c=",ddy
end	

22: begin		;poloidal mesh
    title=ename(choice)
    datax=spdata(*,*,0)
    dataz=spdata(*,*,1)
    multiple,lsp2,lst,datax,dataz,title
end	

23: begin		;b-field
    title=ename(choice)
    x=spdata(*,*,0)
    y=spdata(*,*,1)
    f=spdata(*,*,2)
    conplot,lsp2,lst,x,y,f,title
end	

24: begin		;Jacobian
    title=ename(choice)
    x=spdata(*,*,0)
    y=spdata(*,*,1)
    f=spdata(*,*,3)
    conplot,lsp2,lst,x,y,f,title
end	

25: begin		;icurrent
    title=ename(choice)
    x=spdata(*,*,0)
    y=spdata(*,*,1)
    f=spdata(*,*,4)
    conplot,lsp2,lst,x,y,f,title
end	

26: begin		;zeta2phi
    title=ename(choice)
    x=spdata(*,*,0)
    y=spdata(*,*,1)
    f=spdata(*,*,5)
    conplot,lsp2,lst,x,y,f,title
end	

27: begin		;delb
    title=ename(choice)
    x=spdata(*,*,0)
    y=spdata(*,*,1)
    f=spdata(*,*,6)
    conplot,lsp2,lst,x,y,f,title
end	

28: begin		;b-field (theta)
    title='b-field (theta) at psi=isp'
    x=tdata
    y1=fltarr(lst+1)
    y2=fltarr(lst+1)
    y1(0:lst-1)=spdata(isp,*,2)
    y1(lst)=y1(0)
    y2=y1
    panel1,x,y1,x,y2,title
end	

29: begin		;Jacobian
    title='Jacobian spdata & (gq+I)/B^2 at psi=isp'
    x=tdata
    y1=fltarr(lst+1)
    y2=fltarr(lst+1)
; J data    
    y1(0:lst-1)=3.8*spdata(isp,*,3)
; J=(gq+I)/B^2
    y2(0:lst-1)=(pdata(isp,11)*pdata(isp,10)+spdata(isp,*,4))/(spdata(isp,*,2)*spdata(isp,*,2))
    y1(lst)=y1(0)
    y2(lst)=y2(0)
    panel1,x,y1,x,y2,title
end	

30: begin		;error of spline cos and sin
    title='error of spline cos and sin'
    x=indgen(lsp1)
    y1=pdata(nind0:nind1,28)
    y2=pdata(nind0:nind1,29)
    panel1,x,y1,x,y2,title
end	

31: begin	;plot to PS
    set_plot,'PS'
    print, 'file name=?'
    title='a'
    read,title
    if title ne '' then device, file=title,/color
end	

32: begin       ;change window size
    size_x=1
    size_y=1
    print, 'size_x,size_y=?'
    read,size_x,size_y
    window, plotide, TITLE='equilibrium plot', xsize=size_x,ysize=size_y
end	

33: begin	;change isp
    print, isp
    print, 'new isp=? max=',lsp2-1
    read, isp
end	

34: begin	     ; radial domain
    psi0=0
    psi1=1
    print, 'psin0,psin1=?'
    read,psin0,psin1

    ; sane values
    if psin0 lt 0 then psin0=0
    if psin1 lt 0 then psin1=1
    if psin0 gt 1 then psin0=0
    if psin1 gt 1 then psin1=1

    ; normalize poloidal flux so we use values above
    psin=pdata(*,0)/pdata(lsp1-1,0)

    print,'psi(end)=',pdata(lsp1-1,0)

; find indice range corresponding to user specified domain
    i=0
    while psin(i) le psin0 do begin i=i+1 & endwhile
    nind0=i-1  ; lower boundary index
    while psin(i) lt psin1 do begin i=i+1 & endwhile
    nind1=i    ; upper boundary index

    print,'nind0,nind1=',nind0,nind1
end 

35: begin            ; Quit
    set_plot,'X'
    wdelete                     ; Closes plot windowse
    widget_control, event.top,/destroy
end

endcase
end
;****************************************************************
; panel1 line plot
pro panel1,x1,y1,x2,y2,title
common idle,plotide,ncolor,layer,cvalue,truecolor,pseudo_color,true_color

if !D.name eq 'X' then wset,plotide
!noeras=0
!linetype=0
!p.color=0
!mtitle=title
set_viewport,0.15,0.85,0.15,0.85

xmax=max(x1)
xmin=min(x1)
ymax=max([y1,y2])
ymin=min([y1,y2])
if ymax eq ymin then print,ymax

set_xy,xmin,xmax,ymin,ymax
!linetype=0
plot,x1,y1,charsize=3
!linetype=1
!p.color=truecolor(ncolor*3/4)
oplot,x2,y2

if !D.name eq 'PS' then begin
    device,/close
    set_plot,'X'
    truecolor=true_color
endif

return
end
;****************************************************************
; two-panel plot
pro panel2,x1,y1,x2,y2,title1,title2
common idle,plotide,ncolor,layer,cvalue,truecolor,pseudo_color,true_color

if !D.name eq 'X' then wset,plotide
!noeras=0
!linetype=0
!p.color=0

!mtitle=title1
set_viewport,0.15,0.95,0.58,0.95

print,'x1',x1

xmax=max(x1)
xmin=min(x1)
ymax=max(y1)
ymin=min(y1)

print,'xmax',xmax

if ymax eq ymin then print,ymax
set_xy,xmin,xmax,ymin,ymax
plot,x1,y1,/ystyle,charsize=3

!p.color=truecolor(ncolor*3/4)
!noeras=1
!mtitle=title2
set_viewport,0.15,0.95,0.1,0.47
;set_viewport,0.15,0.95,0.58,0.95
xmax=max(x2)
xmin=min(x2)
ymax=max(y2)
ymin=min(y2)
if ymax eq ymin then print,ymax
set_xy,xmin,xmax,ymin,ymax
plot,x2,y2,/ystyle,charsize=3,psym=2
print,x2,y2

if !D.name eq 'PS' then begin
    device,/close
    set_plot,'X'
    truecolor=true_color
endif

return
end

;****************************************************************
; multiple curves plot
pro multiple,lsp,lst,datax,dataz,title
common idle,plotide,ncolor,layer,cvalue,truecolor,pseudo_color,true_color

 if !D.name eq 'X' then wset,plotide
 !noeras=0
 !linetype=0
 !p.color=0
 !mtitle=title

 lp=lsp-1

 xmax=max(datax(0:lp,*))
 xmin=min(datax(0:lp,*))
 ymax=max(dataz(0:lp,*))
 ymin=min(dataz(0:lp,*))
 print, xmin,xmax,ymin,ymax
 xw=0.2+0.9*(xmax-xmin)/(ymax-ymin)
 if !D.name eq 'PS' then xw=xw*0.75 
 set_viewport,0.2,xw,0.05,0.95

 set_xy,xmin,xmax,ymin,ymax
 x=fltarr(lst+1)
 y=fltarr(lst+1)

 for i=0,lp,1 do begin
     x(0:lst-1)=datax(i,*)
     y(0:lst-1)=dataz(i,*)
     x(lst)=x(0)
     y(lst)=y(0)
     if i eq 0 then plot,x,y,/xstyle,/ystyle,charsize=3
     if i ne 0 then oplot,x,y
 end		
 x=fltarr(lp+1)
 y=fltarr(lp+1)

 for j=0,lst-1,1 do begin
     x=datax(0:lp,j)
     y=dataz(0:lp,j)
     oplot,x,y
 end		

 if !D.name eq 'PS' then begin
     device,/close
     set_plot,'X'
     truecolor=true_color
 endif
 return
end
;****************************************************************
pro conplot,lsp,lst,xx,yy,ff,title

common idle,plotide,ncolor,layer,cvalue,truecolor,pseudo_color,true_color

 if max(ff) eq min(ff) then begin
     print,max(ff)
     return
 endif

 if !D.name eq 'X' then wset,plotide
 !noeras=0
 !linetype=0
 !p.color=0
 !mtitle=title
 !noeras=0
 !x.style=1
 !y.style=1

 x=fltarr(lsp,lst+1)
 y=fltarr(lsp,lst+1)
 f=fltarr(lsp,lst+1)
 x(*,0:lst-1)=xx
 x(*,lst)=xx(*,0)
 y(*,0:lst-1)=yy
 y(*,lst)=yy(*,0)
 f(*,0:lst-1)=ff
 f(*,lst)=ff(*,0)

 xmax=max(x)
 xmin=min(x)
 ymax=max(y)
 ymin=min(y)
 xw=0.2+0.9*(xmax-xmin)/(ymax-ymin)
 set_viewport,0.2,xw,0.05,0.95
 set_xy,xmin,xmax,ymin,ymax

 zmax=max(f)
 zmin=min(f)
 flevel=zmin+(zmax-zmin)*indgen(layer)/(layer-1)
 contour,f,x,y,nlevels=layer,c_colors=cvalue,max_value=zmax,min_value=zmin,levels=flevel,/fill,charsize=3

 if !D.name eq 'PS' then begin
     device,/close
     set_plot,'X'
     truecolor=true_color
 endif
 
 return
end
;****************************************************************
pro eplt

common startup,number_plot,fpath,color_max,color_value,direct_color
common idle,plotide,ncolor,layer,cvalue,truecolor,pseudo_color,true_color,nind0,nind1
common edata,ename,nrplot,nplot,lsp1,lsp2,lst,isp,pdata,tdata,spdata

 plotide=number_plot
 true_color=color_value
 truecolor=true_color
 pseudo_color=direct_color
 ncolor=color_max
 layer=color_max
 cvalue=intarr(layer)
 cvalue=truecolor(indgen(layer)*ncolor/layer)


; 8-bits color
; cvalue=indgen(layer)*ncolor/(layer-1)

; default window
 set_plot,'X'
 !linetype=0
 !p.thick=1
 !p.background=truecolor(ncolor/2)
 !p.color=truecolor(ncolor+1)
 window, plotide, title='equilibrium plot', xsize=900,ysize=900

 openr, plotide, 'equilibrium.out'

; # of 1D radial plots and radial points
 readf, plotide, nrplot,lsp1
 pdata=fltarr(lsp1,nrplot+1)
 readf, plotide, pdata
;0: poloidal flux function
;1: nomalized toroidal flux function
;2: minor radius
;3: major radius
;4: Te
;5: -d(ln(Te))/dr
;6: ne
;7: -d(ln(ne))/dr
;8: Ti
;9: -d(ln(Ti))/dr
;10: ni
;11: -d(ln(ni))/dr
;12: Tf
;13: -d(ln(Tf))/dr
;14: nf
;15: -d(ln(nf))/dr
;16: zeff
;17: toroidal roation
;18: radial electric field
;19: q
;20: d(ln(q))/dr
;21: g
;22: p
;23: r: rpsi
;24: toroidal flux: torpsi
;25: radial grid: rgpsi
;26: inversion of torpsi: psitor
;27: inversion of rgpsi: psirg
;28: error of spline cos
;29: error of spline sin

; # of 2D poloidal plots
 readf, plotide, nplot,lsp2,lst
 spdata=fltarr(lsp2,lst,nplot+2)
 readf, plotide, spdata
; spdata
;0: x
;1: z
;2: b
;3: J
;4: i
;5: zeta2phi
;6: del

; sampling one flux surface using isp
 pi=4.0*atan(1.0)
 tdata=fltarr(lst+1)
 for j=0,lst do begin
     tdata(j)=2.0*pi*j/lst
 endfor
 isp=lsp2-1
 print,'lsp1=',lsp1,'lsp2=',lsp2,'lst=',lst,'isp=',isp

 close, plotide

; initialize plot range
 nind0=0
 nind1=lsp1-1

; make a menu
 ename=strarr(36)
 ename=["minor","major","Te","ne","Ti","ni","Tf","nf","Zeff","rotation","E_r",$
  "q:psi-r","s:psi-1","g:psi-r","p:psi-r","tor:psi-r","r:psi-tor","rg:psi-r","psi(tor)",$
  "psi(rg)","q(psi)","psi-r:rg","spline mesh","b-field","Jacobian","icurrent","zeta2phi",$
  "delb","b(theta)","J(theta)","spline error","PS file","window size","isp","radial dom.","Exit"]

 xmenu,ename,BASE=sbase,SPACE=5,title='equilibrium panel',column=4,xpad=5,ypad=5
 widget_control,sbase,/realize
 xmanager,"equilibrium",sbase,/NO_BLOCK

end
;********************************************************************************
