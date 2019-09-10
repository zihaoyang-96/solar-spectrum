;NAME: INTEGRAL_EMISS_CALC
;PURPOSE: to calculate the density sensitive infrared forbidden line pair ratio considering line-of-sight effect
;KEYWORDS: nwvl1, nwvl2: the wavelength number of the specified line for the emissivity of one particular ion
;                       for Fe XIII 10747, 10798 they are 53416 and 53427
;NOTES: The procedure first calculates the emissivity of the lines at different heights, considering photo-excitation
;   (which is important for forbidden lines), using an electron density model from Gibson et al. (1999) and a temperature
;   model (below 1.05 solar radii logT=6.14, over the height the temperature follows the model from Vásquez et al. (2003)).
;       Then for each height the procedure will integrate the emissivities along the LOS over a distance of ~10 solar radii.
;       The final step is to get the line ratio.

;       The emissivity calculated by CHIANTI is hv N_j A_ji (v: \mu; N_j: the partition of the specified ion at an energy level,
;            i.e. N_j=N(X_j^+m)/N(X^+m); A_ji: Einstein A value), the unit is erg/s; but here we use the keyword /no_de in emiss_calc.pro
;           to set the obtained emissivity with an unit of s^-1.
;
;       From emissivity to intensity: I = integral (emiss * N_e * ioneq * dh), ioneq is the ion percentage from ionization equilibrium. 
;           This will be cancelled when deriving the line ratio, so we don't need to include it.
; 
;AUTHOR: Zihao Yang at Peking University, Sept. 7, 2019.


pro integral_emiss_calc, nwvl1, nwvl2
;nwvl1=53416 & nwvl2=53417

;density model
r=findgen(420)/100+1.0
a=3.6 & b=15.3 & c=0.99 & d=7.34 & e=0.365 & f=4.31
; a=77.1 & b=31.4 & c=0.954 & d=8.3 & e=0.55 & f=4.63
ndens=(a*r^(-b)+c*r^(-d)+e*r^(-f))*1e8

;temperature model
t0=8e5 & a1=0.1 & b1=0.33 & a2=0.55 & b2=6.6
t=t0*(a1+1)/(a1+b1*r^a2+(1-b1)*r^(-b2))
t[0:50]=10^6.14

em1074=dblarr(420)
em1079=dblarr(420)
for ii=0, 419 do begin
    emiss1=emiss_calc('fe_13',temp=alog10(t[ii]),dens=alog10(ndens[ii]),rphot=r[ii],/no_de,radtemp=5900)
    em1074[ii]=emiss1[nwvl1].em
    em1079[ii]=emiss1[nwvl2].em
end


int1074=dblarr(41)
int1079=dblarr(41)
for jj=0,40 do begin
    h=fix(sqrt((r[jj])^2+5^2)*100)/100.
    loc=where(r eq h)
    for kk=jj, loc[0]-1 do begin
        int1074[jj]=int1074[jj]+em1074[kk]*ndens[kk]*0.01*6.957/2.*1e8
        int1079[jj]=int1079[jj]+em1079[kk]*ndens[kk]*0.01*6.957/2.*1e8
    endfor
endfor
int1074=int1074*2
int1079=int1079*2
save, filename='integral_emiss.sav',em1074,em1079,int1074,int1079
end




