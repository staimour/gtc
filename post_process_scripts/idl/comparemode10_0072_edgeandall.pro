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
  plot,xx,yr,charsize=2,/ystyle
  !p.color=truecolor(ncolor*3/4)
  oplot,xx,yi
  
  ; panel 2: mode amplitude history
  ya=sqrt(yr*yr+yi*yi)
  if smoothnum eq 0 then begin
    ypow=alog10(ya)
  endif else begin
    ypow=Smooth(alog10(ya),smoothnum);alog10(ya)
  endelse
  
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
  if smoothnum eq 0 then begin
    yr=yr/exp(gamma*xpow)
    yi=yi/exp(gamma*xpow)
  endif else begin
    yr=Smooth(yr/exp(gamma*xpow),smoothnum);yr/exp(gamma*xpow)
    yi=Smooth(yi/exp(gamma*xpow),smoothnum);yi/exp(gamma*xpow)
  endelse
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
  
  ymax_c=max(yp(where(xp gt 0)));
  xmax_c=xp(where(yp eq ymax_c));
  ymax_c2=max(yp(where(xp lt 0)));
  xmax_c2=xp(where(yp eq ymax_c2));
  print,'frequency=',xmax_c,' amplitude=',ymax_c
  print,'frequency=',xmax_c2,' amplitude=',ymax_c2
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