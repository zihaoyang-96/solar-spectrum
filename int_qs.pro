; NAME: int_qs
; PURPOSE: to calculate the change of line intensity along the solar radial distance for corona above quiet sun (QS) region,
;          using the density model of Gibson et al.(1999) and assuming the temperature below 1.5R_sun is logt=6.2,
;          and after that using the temperature model mentioned in Del Zanna et al.(2018)
;         the EM is assumed at logem=27.8 when the electron density is about logn=9 and changes with the square of electron density.
; PARAMETERS: wmin: the minimum wavelength
;             wmax: the maximum wavelength
;             line: the line name in the form of i.e. 'fe_13'
;             output: the name of the output sav file and txt file, suggesting using the line wavelength, i.e. '10747'
;             num: the atomic number of the element, i.e. 26 (for Fe)
; AUTHOR: Written by Zihao Yang, Peking University, May, 2018


pro int_h, wmin, wmax,line,output,num

r=findgen(15)/10+1
a=3.6
b=15.3
c=0.99
d=7.34
e=0.365
f=4.31
n_e=(a*r^(-b)+c*r^(-d)+e*r^(-f))*1e8
logn=alog10(n_e)

t0=8e5
a1=0.1
b1=0.33
a2=0.55
b2=6.6
t=t0*(a1+1)/(a1+b1*r^a2+(1-b1)*r^(-b2))
t[0:5]=10^6.2

read_abund,'/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2009_asplund.abund',abund,ref
ab=abund(num-1)

for i=0,5 do begin
  ch_synthetic,wmin,wmax,sngl_ion=line,density=n_e[i],/photons,logt_isothermal=6.2,logem_isothermal=(27.8+(logn[i]-logn[0])*2),rphot=r[i],radtemp=6000,/all,ioneq_name=!ioneq_file,save_file='test'
  restgen,file='test.genx',str
  openw,lun,strtrim(output)+'.txt',/get_lun,/append
  printf,lun,((str.int)[0])*ab,r[i]
  free_lun,lun
endfor
for i=6,14 do begin
  ch_synthetic,wmin,wmax,sngl_ion=line,density=n_e[i],/photons,logt_isothermal=alog10(t[i]),logem_isothermal=(27.8+(logn[i]-logn[0])*2),rphot=r[i],radtemp=6000,/all,ioneq_name=!ioneq_file,save_file='test'
  restgen,file='test.genx',str
  openw,lun,strtrim(output)+'.txt',/get_lun,/append
  printf,lun,((str.int)[0])*ab,r[i]
  free_lun,lun
endfor
readcol,strtrim(output)+'.txt',int,r
filesave= strtrim(output)+'.sav'
save,filename=filesave,int,r
end
