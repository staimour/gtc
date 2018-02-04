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