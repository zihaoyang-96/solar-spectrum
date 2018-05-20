; NAME: int_ar
; PURPOSE: to calculate the change of line intensity along the radial distance in active region (AR) corona.
;         and assuming the temperature below 1.4R_sun is logt=6.32, and after that using the temperature model of QS with a factor of 1.7 multplied.
;
; PARAMETERS: wmin: the minimum wavelength
;             wmax: the maximum wavelength
;             line: the line name in the form of i.e. 'fe_13'
;             output: the name of the output sav file and txt file, suggesting using the line wavelength, i.e. '10747'
;             num: the atomic number of the element, i.e. 26 (for Fe)
; AUTHOR: Written by Zihao Yang, Peking University, May, 2018


pro int_ar, wmin, wmax,line,output,num

r=findgen(15)/10+1

n_e=1e8*(2.99*r^(-16)+1.55*r^(-6))*5
logn=alog10(n_e)

t0=8e5
a1=0.
b1=0.47
a2=0.7
b2=6.6
t=t0*(a1+1)/(a1+b1*r^a2+(1-b1)*r^(-b2))*1.7

read_abund,'/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2009_asplund.abund',abund,ref
ab=abund(num-1)

for i=0,4 do begin
  ch_synthetic,wmin,wmax,sngl_ion=line,density=n_e[i],/photons,logt_isothermal=6.32,logem_isothermal=(10.4+(logn[i])*2),rphot=r[i],radtemp=6000,/all,ioneq_name=!ioneq_file,save_file='test'
  restgen,file='test.genx',str
  openw,lun,strtrim(output)+'.txt',/get_lun,/append
  printf,lun,((str.int)[0])*ab,r[i]
  free_lun,lun
endfor
for i=5,14 do begin
  ch_synthetic,wmin,wmax,sngl_ion=line,density=n_e[i],/photons,logt_isothermal=alog10(t[i]),logem_isothermal=(10.4+(logn[i])*2),rphot=r[i],radtemp=6000,/all,ioneq_name=!ioneq_file,save_file='test'
  restgen,file='test.genx',str
  openw,lun,strtrim(output)+'.txt',/get_lun,/append
  printf,lun,((str.int)[0])*ab,r[i]
  free_lun,lun
endfor
readcol,strtrim(output)+'.txt',int,r
filesave= strtrim(output)+'.sav'
save,filename=filesave,int,r
end
