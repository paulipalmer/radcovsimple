import numpy as np
import sys
from subs import *

#---Constants
from datavar import *

#---Read input data
filename    = 'input.file'
nheaderline = 2
dval = readinput(filename,nheaderline)

#---Print some output to screen
writetoscreen(dval,'zero')

#print('CLOUD PROPERTIES')
#print('LOW | MIDDLE | HIGH')
#print('pc',dval.pcn)
#print('dpc',dval.dpcn)
#print('clt',dval.cltn)
#print('cla',dval.clan)
#print('ec',dval.ecn)
#print('fr',dval.frn)
#sys.exit()   


#---Start calculation
# Mean daily incoming solar flux at the TOA
flux = dval.amu * dval.day * dval.s0

tg   = dval.t0

# pressure, temperature, altitude, numbers of levels, level of tropopause
p, t, z, it = pgrid(dval)

kk  = 0
er  = 1
tgb = 0
gb  = 0

if dval.itf == 0: const_nmax = 0

while kk <= const_nmax and er > const_ep:
    kk = kk + 1
    print('ITERATION = ', kk)
    ales = 58.1717-6938.67/tg-5.5189*np.log(tg)
    es   = np.exp(ales)
    beta = 0.634*dval.rho*es/dval.wh2o - 1.
    print('tg = ', tg, ' beta = ', beta)

    # Set up pressure grid
    tt = zet(dval,p,t,z,it)



    sum1 = 0. # net incoming solar radiation at the TOA
    sum2 = 0. # Net outgoing terrestrial flux at the TOA
    sum3 = 0. # 
    sum4 = 0. # Net downward solar flux at ground
    sum5 = 0. # Net upward terrestrial flux at ground

    # Loop over cloud levels: low, middle, high cloud
    for ll in np.arange(len(dval.pcn)):

        if ll == 0: print('LOW CLOUD')
        elif ll == 1: print('MIDDLE CLOUD')
        elif ll == 2: print('HIGH CLOUD')

        icb = dval.icbn[ll]
        ic  = dval.icn[ll]
        pc  = dval.pcn[ll]*dval.pg
        dpc = dval.dpcn[ll]
        clt = dval.cltn[ll]
        cla = dval.clan[ll]
        ec  = dval.ecn[ll]

        # Computes net solar flux at the top and base of the atmosphere
        flin, fsolg = solar(dval,pc,clt,cla,flux)

        # Computes the net terrestrial flux at the top and base of the atmosphere
        flout, firg, tbs, tbc, ebs, ebc = terre(dval,tt,t,p,z,ic,icb,beta,ec)

        sum1 = sum1+dval.frn[ll]*flin
        sum2 = sum2+dval.frn[ll]*flout
        sum4 = sum4+dval.frn[ll]*fsolg
        sum5 = sum5+dval.frn[ll]*firg
        
        
    ga  = sum2 - sum1
    dt  = -ga*(tg-tgb)/(ga-gb)
    tgb = tg

    
    if kk ==1: dval.t0 = tgb + 10.
    if kk > 1: dval.t0 = tgb + dt

    tg = dval.t0

    gb = ga

    if dval.itf == 0.: tg = dval.t0

    er = abs((sum2-sum1)/sum1)

    print('KK SUM1 SUM2 DT TG ER')
    print(kk, sum1, sum2, dt, tg, er)

    print('NET IR OUT = ', sum2, ' NET SOL IN = ', sum1)
    print('NET SOLG DOWN = ', sum4, ' NET IRG UP = ', sum5)

    
    egst4 = dval.eg*const_sb*tg**4.

    print('GROUND EMISSION = ', egst4)
    if dval.itf == 0: print('ITF = 0: SURFACE TEMPERATURE SPECIFIED')
    
    if er > const_ep and dval.itf == 1: print('No convergence...')




    
alti(dval,icb,ic,p,z,t,beta)
        
        
# Note to self
# archive meta data for run in output file
