;*******************************************************************************
pro GTC_event,event

common startup,number_plot,fpath,color_max,color_value,direct_color

widget_control,event.id,get_uvalue=choice

case choice of

0: begin            
        widget_control, event.top,/destroy
	exit
end 

1: begin            
        number_plot=number_plot+1
        history
end 

2: begin            
        number_plot=number_plot+1
        snapshot
end 

3: begin            
        number_plot=number_plot+1
        eplt
end 

4: begin            
        number_plot=number_plot+1
        rtime
end 

5: begin            ;movie
        number_plot=number_plot+1
        envelope
end 

6: begin            ;gauss
        number_plot=number_plot+1
        spectrum
end 

7: begin            ;old snap
        number_plot=number_plot+1
        splt
end 

8: begin            ;old history
        number_plot=number_plot+1
        hplt
end 

9: begin            ;colormap

        set_plot,'X'
	!linetype=0
	!p.background=truecolor(color_max+1)
	!p.color=truecolor(color_max)

	x=[0.0,0.0]
	y=[0.0,1.0]
        xmin=0.0
        xmax=float(color_max+2)
        ymin=0.0
        ymax=1.0
        set_xy,xmin,xmax,ymin,ymax
	plot,x,y

;        !p.thick=10

	for i=0,color_max+1 do begin
		xa=float(i)
		x=[xa,xa]
		y=[0.0,1.0]
		!p.color=color_value(i)
		oplot,x,y
	endfor

end 

endcase
end

;*******************************************************************************
common startup,number_plot,fpath,color_max,color_value,direct_color

; read color map
	openr,1,'${HOME}/idl/color.dat'	
	color_max=1
	readf,1,color_max
	red=intarr(color_max+2)
	green=intarr(color_max+2)
	blue=intarr(color_max+2)
	readf,1,red,green,blue
	close,1	
	
	color_value=indgen(color_max+2)
	direct_color=intarr(color_max+2,3)
	direct_color[*,0]=red
	direct_color[*,1]=green
	direct_color[*,2]=blue

nbyte=3
;print, '# of color bytes'
;read, nbyte
;print,nbyte	

;true color
	if nbyte eq 3 then color_value=red+256l*(green+256l*blue)

; load 8-bits color table
	;if nbyte eq 1 then tvlct,red,green,blue
	if nbyte eq 1 then tvlct,direct_color

!p.thick=1


Number_plot=0
fpath='.'

pname=strarr(10)
pname=["Exit idl","History","Snapshot","Equilibrium","Radial-time","old snap","old history","colormap"]

xmenu,pname,Base=pbase,Space=10,Title='GTC_IDL',xpad=20,ypad=20


Widget_Control,pbase,SCR_Xsize=160,Default_font='times-roman',/REALIZE
xmanager,"GTC",pbase,/NO_BLOCK
end
