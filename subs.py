import numpy as np
import sys

#---Constants
from datavar import *

# Simple class to contain data
class datavar:
    def __init__(self,s0,amu,\
                 day,rg,eg,pg,\
                 glapz,pt,cc,\
                 wco2,wh2o,wo3,wch4,rho,\
                 pcn, dpcn,cltn,clan,ecn,frn,\
                 t0,itf):
        self.s0    = s0
        self.amu   = amu
        self.day   = day
        self.rg    = rg
        self.eg    = eg
        self.pg    = pg
        self.glapz = glapz
        self.pt    = pt
        self.cc    = cc
        self.wco2  = wco2
        self.wh2o  = wh2o
        self.wo3   = wo3
        self.wch4  = wch4
        self.rho   = rho
        self.pcn   = pcn
        self.dpcn  = dpcn
        self.cltn  = cltn
        self.clan  = clan
        self.ecn   = ecn
        self.frn   = frn
        self.t0    = t0
        self.itf   = itf
        self.icn   = np.zeros(3)
        self.icbn  = np.zeros(3)
        

def readinput(filename,nheaderline):
    f = open(filename,'r')
    # Read informative header
    for ii in np.arange(nheaderline):
        inline = f.readline()
        #print(inline.strip())

    vval  = []

    # Read input
    while 1:
        inline = f.readline()
        if not inline: break
        cols = inline.rsplit('=')
        vval.append(np.float(cols[1]))
    f.close()

    dataval = datavar(vval[0],vval[1],\
                      vval[2],vval[3],vval[4],vval[5],\
                      vval[6],vval[7],vval[8],\
                      vval[9],vval[10],vval[11],vval[12],vval[13],\
                      vval[14:16+1],vval[17:19+1],vval[20:22+1],\
                      vval[23:25+1],vval[26:28+1],vval[29:31+1],\
                      vval[32],vval[33])

    return dataval



def writetoscreen(dataval,label):
    if label=='zero':
        print('--SIMPLE CLIMATE MODEL---')
    elif label == 'first':
        keys = list(dataval.__dict__.keys())
        vals = list(dataval.__dict__.values())
        for ii in np.arange(len(keys)): print(keys[ii],'==',vals[ii])


def pgrid(dval):
    # Sets up pressure grid and locates the cloud and tropopause levels.

    p = np.zeros(imax)
    t = np.zeros(imax)
    z = np.zeros(imax)
    
    # Pressure grid interval set equal to 10 hPa
    dp = dval.pg/np.float(imax-1)

    for ii in np.arange(imax):

        p[ii] = dval.pg - dp * np.float(ii)

        for jj in np.arange(len(dval.pcn)):
            # Pressure at top of cloud type N
            pc  = dval.pcn[jj]*dval.pg
            # Presure at bottom of cloud type N
            pcb = pc + dval.dpcn[jj]

            # icn corresponds to level of top of cloud type N
            if p[ii] == pc:  dval.icn[jj]  = ii
            # icn corresponds to level of bottom of cloud type N
            if p[ii] == pcb: dval.icbn[jj] = ii

        # it corresponds to level of tropopause
        if p[ii] == dval.pt: it = ii

    return p, t, z, it


def zet(dval,p,t,z,it):
    # Sets up the temperature versus presure profile and the altitude
    # versus pressure profile.

    tg = dval.t0
    
    alfa = dval.glapz*const_r/(const_g*const_amw)

    zt = tg*(1.-(dval.pt/dval.pg)**alfa)/dval.glapz

    imax1 = imax - 1

    for ii in np.arange(imax1):
        pi = p[ii]
        pr = pi/dval.pg
        if pi >= dval.pt: z[ii] = tg*(1.-pr**alfa)/dval.glapz
        t[ii] = tg-dval.glapz*z[ii]
        if pi < dval.pt:
            t[ii] = tg-dval.glapz*zt
            z[ii] = zt-np.log(pi/dval.pt)*alfa*t[ii]/dval.glapz
        if z[ii] == zt: it = ii

    # tropopause temperature
    tt      = t[it]
    
    z[imax-1] = 100.
    t[imax-1] = tg - dval.glapz*zt

    return tt


def solar(dval,pc,clt,cla,flux):
    # Computes net solar flux at the top and base of the atmosphere.

    pcg = pc/dval.pg

    r3g = data_r3 * (1. - data_gs)*clt

    rc  = r3g/(2.+r3g)

    ctr = 1.-rc-cla

    twp = 1.-0.11*(dval.wh2o/dval.amu)**0.31
    twm = 1.-0.11*(1.66*dval.wh2o)**0.31

    tvp = 1.-0.021*dval.wo3/dval.amu
    tvm = 1.-0.021*1.66*dval.wo3
    tuv = 1.-0.023*(dval.wo3/dval.amu)**0.18

    tcpp = 1.-0.015*(dval.wco2*pcg/dval.amu)**0.263

    tcpm = 1.-0.015*(dval.wco2*pcg*1.66)**0.263
    tcp  = 1.-0.015*(dval.wco2/dval.amu)**0.263
    tcm  = 1.-0.015*(1.66*dval.wco2)**0.263
    tcpd = 1.-0.015*(dval.wco2*(1.-pcg)*1.66)**0.263

    ap = tuv*tvp*tcpp
    am = tvm*tcpm

    bp = data_tr*twp*tcp*tvp*tuv
    bm = data_tr*twm*tcm*tvm

    aa = ap*am
    bb = bp*bm

    bc = twm*tcpd
    bd = aa*ctr*ctr*bc*bc

    ac = 1.-aa*rc - bd*dval.rg
    fsc = flux * ac

    vas = 1.-(1.-data_tr)*tuv*tvp*tvm-bb*dval.rg

    fsnc = flux*vas

    flin = dval.cc*fsc+(1.-dval.cc)*fsnc

    tstars = bp*(1.-dval.rg)
    tstarc = ap*ctr*bc*(1.-dval.rg)

    fsolg = (1.-dval.cc)*flux*tstars+dval.cc*flux*tstarc
    
    return flin, fsolg


def alti(dval,icb,ic,p,z,t,beta):

    tg = dval.t0

    f = open('profile.out', 'w')
    f.write('VERTICAL STRUCTURE\n')
    f.write('Z    P    T   TAP   TACB    TAC    TAUD\n')
    
    for ii in np.arange(imax):
        tap  = taus(0,ii,dval,p,z,t,beta)
        tacb = taus(ii,icb,dval,p,z,t,beta)
        tac  = taus(ii,ic,dval,p,z,t,beta)
        taud = taus(ii,imax-1,dval,p,z,t,beta)

        hp = htau(t[ii],tap,tg)
        hm = htau(t[ii],taud,tg)

        ft = -0.0045*t[1]+2.3
        
        strout = "%8.2f %8.2f %8.2f %8.4f %8.4f %8.4f %8.4f \n" % \
                 (z[ii], p[ii], t[ii], tap, tacb, tac, taud)
        f.write(strout)
        
    f.close()
        

def htau(uset,usetau,tg):
    # Computes integrand H(r)
    return np.exp(-usetau)*((uset/tg)**4.)


def oz(z):
    # Computer Green's function G(z)
    an  = 1. + np.exp(-data_b/data_c)
    den = 1. + np.exp((z-data_b)/data_c)
    return an/den
                  

def taus(i1,i2,dval,p,z,t,beta):
    p1 = p[i1]
    p2 = p[i2]
    pg = p[0]
    z1 = z[i1]
    z2 = z[i2]

    ozf1 = oz(z1)
    ozf2 = oz(z2)

    ft = -0.0045*t[0]+2.3

    pr  = abs((p1-p2)/pg)
    pr5 = (p1/pg)**(1.+beta)-(p2/pg)**(1.+beta)
    pr5 = abs(pr5)

    wh  = dval.wh2o*pr5
    wc  = dval.wco2*pr
    wo  = 0.00214*dval.wo3*abs(ozf1-ozf2)
    wch = dval.wch4*pr

    tah  = 0.63*wh**0.52
    tac  = 0.14*wc**0.22
    tao  = 2.51*wo**0.62
    tach = 2.51*wch**0.75

    return 1.66*ft*(tah+tac+tao+tach)



def hint(ip,i1,im,dval,p,z,t,beta,tg):
    # Evaluates integrals using the trapezoidal rule
    
    taup = np.zeros(imax)
    
    for ii in np.arange(imax):
        tau = taus(ip,ii,dval,p,z,t,beta)
        taup[ii] = tau

    i2  = i1 + 1
    im1 = im - 1

    dp = taup[i2]-taup[i1]
    tp = t[i1]
    taupp = taup[i1]
    hp = htau(tp,taupp,tg)
    sum1 = -dp*hp/2.

    # Note +1 is to ensure that the counter gets to im1
    for ii in np.arange(i2,im1+1):
        dp = taup[ii+1] - taup[ii-1]
        tp = t[ii]
        taupp = taup[ii]
        hp = htau(tp,taupp,tg)
        sum1 = sum1 - dp*hp/2.

    dp    = taup[im]-taup[im1]
    tp    = t[im]
    taupp = taup[im]
    hp    = htau(tp,taupp,tg)

    fip = sum1 - dp*hp/2.
    print('fip = ', fip)
    return fip

    
def terre(dval,tt,t,p,z,ic,icb,beta,ec):
    # Computes the net terrestrial flux at the top and base of the atmosphere

    tg = dval.t0
    
    tc = t[ic]

    taug  = taus(0,imax-1,dval,p,z,t,beta)
    tauc  = taus(ic,imax-1,dval,p,z,t,beta)
    taucb = taus(0,icb,dval,p,z,t,beta)
    
    trg  = np.exp(-taug)
    trc  = np.exp(-tauc)
    trcb = np.exp(-taucb)

    t4  = tg**4.
    tc4 = tc**4.

    print('TG TC TT TAUG TAUC')
    print(tg,tc,tt,taug,tauc)

    fsout = hint(imax-1,1-1,imax-1,dval,p,z,t,beta,tg)
    fcout = hint(imax-1,ic,imax-1,dval,p,z,t,beta,tg)
    fs0   = hint(1-1,1-1,imax-1,dval,p,z,t,beta,tg)
    fcub  = hint(icb,1-1,icb,dval,p,z,t,beta,tg)
    fcdb  = hint(1-1,1-1,icb,dval,p,z,t,beta,tg)
    fcda  = hint(ic,ic,imax-1,dval,p,z,t,beta,tg)

    egst4 = dval.eg*const_sb*t4

    tbs = trg+fsout/dval.eg - (1.-dval.eg)*trg*fs0/dval.eg
    fcl = tbs*egst4
    tbc = trc+(t4/tc4)*fcout/ec+(1.-ec)*(fcub+dval.eg*trcb)*trc*t4/(ec*tc4)
    fac = tbc*ec*const_sb*tc4

    flout = (1.-dval.cc)*fcl+dval.cc*fac

    print('TBS TBC')
    print(tbs,tbc)
    print('IROUT CLEAR = ', fcl,' IROUT CLOUD = ', fac)
    print('IROUT TOTAL = ', flout)

    fbsd = egst4*(-fs0)
    fbsu = egst4+(1.-dval.eg)*fbsd/dval.eg
    fbs  = fbsu-fbsd

    fbcd = egst4*(-fcdb+ec*(tc4/t4)*trcb - (1.-ec)*fcda*trcb)
    fbcu = egst4+(1.-dval.eg)*fbcd/dval.eg
    fbc  = fbcu-fbcd
    ebc  = fbc/(const_sb*t4)
    ebs  = fbs/(const_sb*t4)

    print('EBS EBC')
    print(ebs,ebc)

    firg = (1.-dval.cc)*fbs+dval.cc*fbc

    print('FBS FBC CC FIRG')
    print(fbs,fbc,dval.cc,firg)

    return flout,firg,tbs,tbc,ebs,ebc
    
