;#NAME: EVE_MEDIAN_SPECTRA

;#AUTHOR:Zihao Yang, Peking University, April 1, 2018

;#PURPOSE: Derive the median spectrum of 2 hr EVE data
;           

;#CAUTION: You need to specify the directory where your data are stored, if you want to merge data for longer durations, 
;           you should change the value '720'  of 'a=fltarr(720)' and reset the steps of the loops.

;#OUTPUT: A save file with parameters WVL (the eve wavelength, unchanged) and 
;          IRR (median of the 720 10s observed irradiance for all 5200 wvl bin)

pro eve_median_spectra
infn=file_search("/Research_Data/Schonfeld_2017/2011195/EVS_L2_2011195_*_006_02.fit.gz",count=count)
l2s1=eve_read_whole_fits(infn(0),verbose=verbose)
l2s2=eve_read_whole_fits(infn(1),verbose=verbose)
irr=fltarr(5200)
for m=0,5199 do begin
  a=fltarr(720)
  for i=0,359 do begin
    a[i]=(l2s1.spectrum.irradiance[m])[i]
    endfor
  for j=360,719 do begin
    a[j]=(l2s2.spectrum.irradiance[m])[j-360]
    endfor
  irr[m]=median(a)
  endfor
wvl=l2s1.spectrummeta.wavelength
save,irr,wvl,filename='EVS_2011195_median.sav'
end
