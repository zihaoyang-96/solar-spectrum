; #NAME: SGF_EVE_SPEC

; #AUTHOR:Zihao Yang, Peking University, April 1, 2018

; #PURPOSE: Derive the centroid wavelength, peak intenisty and the uncertainty of measurements of certain line from SDO/EVE spectrum
;           Peak Intensity: erg/s/cm^-2
;           Centroid Wavelength: Angstrom

; #Parameters: wav_min: the minimum wavelength limit
;              wav_max: the maximum wavelength limit
;              lineID: the line name in the form of e.g.' Fe XI'

pro sgf_eve_spec_1, file=file, wav_min, wav_max, lineID=lineID

restore, file=file, /ver

w1=eve_get_wave_bin(wav_min, /l2)
w2=eve_get_wave_bin(wav_max, /l2)

wav=wvl[w1:w2]
irr=irr[w1:w2]

i=w2-w1+1
err=dblarr(i)
;int=dblarr(i)

for j=0,i-1 do begin
  err[j]=irr[j]*0.2
end


res=mpfitpeak(wav, irr, a, nterms=4, /double, /positive)
wave0=a[1]
print,'peak intensity:',a[0]*sqrt(3.14)*a[2]/(6.78e-5)
print,'centroid wavelength:',a[1]

dlambda=eve_get_wave_bin(/inv,/l2,w1+1)-eve_get_wave_bin(/inv,/l2,w1)

loadct,0  & tvlct,255L,0L,0L,2L  & tvlct,0L,255L,0L,3L   & tvlct,0L,0L,255L,4L
window, 1, xs =800, ys =600, xpos=800, ypos=500

plot_io,wav,irr,xrange=[min(wav)-dlambda/2, max(wav)+dlambda/2],yrang=[min(irr)*0.99,max(irr)*1.1]
x=wav
err_plot,x,irr-err,irr+err,width=0.005,thick=1
oplot,x,spline(x,res,x),color=3L

k=indgen(i)
chisq = (1./(k - 4)) * total(((irr[k] - res[k])/err[k])^2)

err_ave=average(err)
print,'err:',err_ave*sqrt(3.14)*a[2]/(6.78e-5)

filesave= 'line'+lineID+'.sav'
save,filename=filesave,a,chisq,err,wave0,res,err_ave
end
