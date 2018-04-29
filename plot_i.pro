pro plot_i,file,line,wvl

restore,file

ymin=min(alog10(int))
ymax=max(alog10(int))

set_plot,'ps'
Device,FileName=strtrim(wvl)+'.eps',XSize=25,YSize=20,/Color,Bits=8,/Encapsul
plot,r,alog10(int),xr=[1,2],yr=[ymin-0.5,ymax+0.5],xstyle=1,xtitle='Radial distance in solar radii',$
  ytitle='Line intensity (log photons !N cm!E-2!N !N sr!E-1!N !N s!E-1!N)',title=strtrim(line)+' '+strtrim(wvl)+' '+'line (photospheric abundance)',$
  charsize=2.5,font=1,thick=3
Device, /close

set_plot,'ps'
Device,FileName=strtrim(wvl)+'_ext.eps',XSize=25,YSize=20,/Color,Bits=8,/Encapsul
plot,r,alog10(int),xr=[1,5],yr=[ymin-0.5,ymax+0.5],xstyle=1,xtitle='Radial distance in solar radii',$
  ytitle='Line intensity (log photons !N cm!E-2!N !N sr!E-1!N !N s!E-1!N)',title=strtrim(line)+' '+strtrim(wvl)+' line (photospheric abundance)',$
  charsize=2.5,font=1,thick=3
Device, /close
end
