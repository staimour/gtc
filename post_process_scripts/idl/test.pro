restore, '/global/u1/c/clau/m=10_all_00725_modehist.sav'
x1=modehist(0:16000,0:1,0,0)
restore, '/global/u1/c/clau/m=10_edge_00725_modehist.sav'
x2=modehist(0:16000,0:1,0,0)

y=[[x1],[x2]]
; panel 1: mode history of real and imaginary components

yrall=y(0:16000,0)
yiall=y(0:16000,1)
yredge=y(0:16000,2)
yiedge=y(0:16000,3)

xtime=indgen(16001)
xx=xtime(0:16000)

xmax=max(xx)
xmin=min(xx)
ymax=max([yrall,yiall,yredge,yiedge])
ymin=min([yrall,yiall,yredge,yiedge])
set_xy,xmin,xmax,ymin,ymax
set_viewport,0.14,0.49,0.6,0.95
!linetype=0
!p.color=truecolor(ncolor+1)
plot,xx,yrall,charsize=2,/ystyle
!p.color=truecolor(ncolor*3/4)
oplot,xx,yiall
oplot,xx,yredge
oplot,xx,yiedge