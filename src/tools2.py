#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 17:37:24 2019

collection of functions dealing with mixing

Z17 => Zhang et al. (2017)
A18 => Armstrong et al. (2019) OR Armstrong 2018, thesis
O16 => O'Neill et al. (2006)
@author: jiedeng
"""
import uncertainties as uct
from uncertainties import umath  
from src.tools import cal_PV
import src.tools as tl
import numpy as np
import scipy.optimize as opt
from scipy.interpolate import interp1d 
import pandas as pd
from src.const import (par_12p5_re, par_12p5_ox, par_25_re, par_25_ox,
                      z17, df_z17, A18_FeO_re, A18_FeO_ox, A3_to_cm3, cm3_to_A3, R)

xlsx = 'db/prev_sum.xlsx'
W = pd.read_excel(xlsx,sheet_name='molar',nrows=4,usecols=range(25),index_col=0)
planets = pd.read_excel(xlsx,sheet_name='planets',nrows=3,skiprows=None,usecols = range(10),index_col=0)

def unpack_par(par,flag):
    """
    unpack the par,

    """
    if flag:
        a = uct.ufloat(par['a'].values, par['aerr'].values)
        b = uct.ufloat(par['b'].values, par['berr'].values)
        W_Fe  = uct.ufloat(par['W_Fe'].values, par['W_Feerr'].values)
        W_Mg  = uct.ufloat(par['W_Mg'].values, par['W_Mgerr'].values)
        W_Si  = uct.ufloat(par['W_Si'].values, par['W_Sierr'].values)
        W_Al  = uct.ufloat(par['W_Al'].values, par['W_Alerr'].values)
        W_Ca  = uct.ufloat(par['W_Ca'].values, par['W_Caerr'].values)
        W_Na  = uct.ufloat(par['W_Na'].values, par['W_Naerr'].values)
        W_K   = uct.ufloat(par['W_K'].values, par['W_Kerr'].values)
        W_Ph  = uct.ufloat(par['W_Ph'].values, par['W_Pherr'].values)
        W_Ti  = uct.ufloat(par['W_Ti'].values, par['W_Tierr'].values)
    else:
#        print("par is",par)
        a = par['a']
        b = par['b']
        W_Fe  = par['W_Fe']
        W_Mg  = par['W_Mg']
        W_Si  = par['W_Si']
        W_Al  = par['W_Al']
        W_Ca  = par['W_Ca']
        W_Na  = par['W_Na']
        W_K   = par['W_K']
        W_Ph  = par['W_Ph']
        W_Ti  = par['W_Ti']
    return a,b,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti


def fo2_cal_single(P,T,PV,r,method='earth',flag=False,fit='fit3'):
    """
    calculate fo2 for a single P,T point
    P : pressure, GPa, P is actually not used, it is included in PV integration info
    T : temperature, K
    PV : term int\{P0 to P0, deltaV}, it requires to get the dV at every single pressure
         and temperature, and then integrate. This could be very time demanding.
    r : Fe3/Fe2 ratio
    method : earth, mars, moon, specify the composition, default is earth
    flag : boolean, uncertainty flag, default is False
    fit3 : W model, default is fit3
    
    Note
    ----

    
    Term int\{P0 to P0, deltaV}, deltaV is at reference temperature, not along geotherm.
    That means, to calcualte this term, for every point in geotherm, P,T, one needs to
    calculate the dV from P0 to P along T, and then integrate.
    Note this is different from just calculate dV along the goetherm, then integrate.
    
    Typically, the 2nd method will underestimate the PV term since goetherm is T0 to T,
    whereas here we need keep T as T.
 
     Nevertheless,  `cal_PV_geo` function in the package solve this problem. The obtained
     PV term along the geotherm is then used as input and kept in the excel eos.xlsx
    """
    
    Fe2 = (1-r)*planets.loc[method]['FeO']
    Fe3 = r*planets.loc[method]['FeO']
    Si = planets.loc[method]['SiO2'];
    Mg = planets.loc[method]['MgO'];
    Al = planets.loc[method]['Al2O3'];
    Ca = planets.loc[method]['CaO'];
    K  = planets.loc[method]['K2O'];
    Na = planets.loc[method]['Na2O'];
    Ph = planets.loc[method]['P2O5'];
    Ti = planets.loc[method]['TiO2'];
    
    a,b,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti = unpack_par(W.loc[fit],flag)
    
    Gterm = np.array(tl.Gr_janaf(T,flag))/R/T
    Wterm = (W_Fe*(Fe2 - Fe3) + W_Mg*Mg + W_Al*Al + W_Si*Si +
             W_Ca*Ca +W_Na*Na + W_K*K + W_Ph*Ph + W_Ti*Ti)/R/T
    tmp   = (PV + Gterm + Wterm)*4
    lnXfe3_Xfe2 = np.log(Fe3/Fe2)*4
    fg = umath.exp(tmp+lnXfe3_Xfe2)
    log10fg = (umath.log10(fg))
    return log10fg

def fo2_cal_vector(P,T,dV,r,method='earth',flag=False,fit='fit3',PV_cal = 'old'):
    """
    calculate XFe3/Xtot as a function at fixed fo2, as a function of P
    -----
    input  list
    fo2 -> absolute fo2
    r   -> (Fe3/XFe)
    """
    PV = np.zeros(P.shape)
    PV_table = pd.read_excel(tl.xlsx,sheet_name='PV',index_col=None)
    if PV_cal == 'old':
        for i in range(len(P)):
            PV[i]  = cal_PV(dV[range(i+1)],P[range(i+1)],T[i],P[i],minP=0)+0.06 # 0.06 is based on correction
    else:
        for i in range(len(P)):
            f = interp1d(PV_table['P(GPa)'],PV_table[PV_cal])
            PV[i]  = f(P[i])
    f_v = np.vectorize(fo2_cal_single,excluded=[3,4,5,6,7])
    return f_v(PV,T,PV,r,method=method,flag=flag,fit=fit)

def X_cal_single(P,T,PV,fo2,method='earth',flag=False,fit='fit3'):
    """
    calculate XFe3/Xtot as a function at fixed fo2, as a function of P
    -----
    input  list
    fo2 -> absolute fo2
    r   -> (Fe3/XFe)
    """
    func = lambda r: fo2_cal_single(P,T,PV,r,method,flag,fit=fit) - fo2
    try:
 #       out = opt.brentq(func, 0.0005, 0.1) # for iw
        out = opt.brentq(func, 0.0005, 0.9) # for iw
        return out
    except:
        print("r =  0.0005, func is", fo2_cal_single(P,T,PV,0.0005,method=method,flag=flag,fit=fit) - fo2)
        print("r =  .99, func is", fo2_cal_single(P,T,PV,.99,method=method,flag=flag,fit=fit) - fo2)

def X_cal_vector(P,T,dV,fo2,method='z17',flag=False,fit='fit3',PV_cal = 'old'):
    """
    calculate XFe3/Xtot as a function at fixed fo2, as a function of P
    -----
    input  list
    fo2 -> absolute fo2
    r   -> (Fe3/XFe)
    
    Note
    -----
    see a detailed explanation of PV in `fo2_cal_single`. Basically, `cal_PV` assume
    PV is along isotherm. Note the input `dV` is usually calcualted with `cal_dV_this_study`
    
    """
    PV = np.zeros(P.shape)
    PV_table = pd.read_excel(tl.xlsx,sheet_name='PV',index_col=None)
    if PV_cal == 'old':
        for i in range(len(P)):
            PV[i]  = cal_PV(dV[range(i+1)],P[range(i+1)],T[i],P[i],minP=0)+0.06 # 0.06 is based on correction
    else:
        for i in range(len(P)):
            f = interp1d(PV_table['P(GPa)'],PV_table[PV_cal])
            PV[i]  = f(P[i])
    f_v = np.vectorize(X_cal_single,excluded=[4,5,6,7,8])
    return f_v(P,T,PV,fo2,method=method,flag=flag,fit=fit)

def cal_PV_geo(P,T,name='Fe_25'):
    """
    calcualte PV for P,T along geotherm
    added by Jie on Oct 23, 2019 
    """
    PV_geo = []
    for i in range(len(P)):
#        print("P is",P[i])
        try:
            Parray = np.linspace(1e-4,P[i])
            Tarray = np.ones_like(Parray)*T[i]
            _, _,_, _, dv = tl.cal_dV_this_study(Parray,Tarray,name=name,flag=False)
        except:
            try:
                Parray = np.linspace(1,P[i])
                Tarray = np.ones_like(Parray)*T[i]
                _, _,_, _, dv = tl.cal_dV_this_study(Parray,Tarray,name=name,flag=False)                   
            except:
                Parray = np.linspace(2,P[i])
                Tarray = np.ones_like(Parray)*T[i]
                _, _,_, _, dv = tl.cal_dV_this_study(Parray,Tarray,name=name,flag=False)    
        PV_geo.append(tl.cal_PV(dv,Parray,T[i],maxP=P[i]))
    return np.array(PV_geo)
        

######################################################
####################### Z17 ###########################

def fo2_Z17(P,T,r,dVdT = 2.92,kmodel = 'K0_old'):
    """
    calculate XFe3/Xtot as a function at fixed fo2, as a function of P
    r is Fe3+/tot Fe
    -----
    input  list
    P,T shoud be scalar 
    """
    x  = r/(1-r)
    tmp = df_z17.loc[dVdT]
    res = np.log(x) - tmp['a'] - tmp['b']/R/T + \
          (20170 + 4.54*(T-1673))*16.6/3*((1+0.241*P)**0.75-1)/R/T - \
          (tmp['c']+dVdT*(T-1673))*tmp[kmodel]/3*((1+4/tmp[kmodel]*P)**0.75-1)/R/T
          #(tmp['c']+dVdT*(t-1673))*36.61/3*((1+4/36.61*25)**0.75-1)/R/t
    fo2 = np.exp(4*res)    
    return np.log10(fo2)


def X_Z17(P,T,fo2,dVdT = 2.92,kmodel = 'K0_old'):
    tmp = df_z17.loc[dVdT]
    res   = 1/4*np.log(10**fo2)+tmp['a'] + tmp['b']/R/T - \
          (20170 + 4.54*(T-1673))*16.6/3*((1+0.241*P)**0.75-1)/R/T + \
          (tmp['c']+dVdT*(T-1673))*tmp[kmodel]/3*((1+4/tmp[kmodel]*P)**0.75-1)/R/T
    x =   np.exp(res)
    r = x/(1 + x)
    return r

######################################################
####################### O16 ###########################

def fo2_O06(P,T,r,method='earth'):
    
    Fe2 = (1-r)*planets.loc[method]['FeO']
    Fe3 = r*planets.loc[method]['FeO']    
    Si = planets.loc[method]['SiO2'];
    Mg = planets.loc[method]['MgO'];
    Al = planets.loc[method]['Al2O3'];
    Ca = planets.loc[method]['CaO'];
    K  = planets.loc[method]['K2O'];
    Na = planets.loc[method]['Na2O'];
    Ph = planets.loc[method]['P2O5'];    
    Ti = planets.loc[method]['TiO2'];

    x = r/(1-r)
    log10fo2 = 4*np.log10(x)-28144/T + 13.95 + 3905*Mg/T - 13359*Ca/T - 14858*Na/T - \
         9805*K/T + 10906*Al/T + 110971*Ph/T - \
         11952*(Fe2 -Fe3)/T + (33122/T - 5.24)*((1+0.241*P)**(3/4)-1) - \
         (39156/T - 6.17)*((1+0.132*P)**(3/4)-1)
    return log10fo2

def X_O06_single(P,T,fo2,method='earth'):
    func = lambda r: fo2_O06(P,T,r,method) - fo2
    try:
 #       out = opt.brentq(func, 0.0005, 0.1) # for iw
        out = opt.brentq(func, 0.0005, 0.9) # for iw
        return out
    except:
        print("r =  0.0005, func is", fo2_O06(P,T,0.0005,method) - fo2)
        print("r =  .99, func is", fo2_O06(P,T,.99,method) - fo2)    

def X_O06(P,T,fo2,method='z17'):
    """
    calculate XFe3/Xtot as a function at fixed fo2, as a function of P
    -----
    input  list
    fo2 -> absolute fo2
    r   -> (Fe3/XFe)
    """

    f_v = np.vectorize(X_O06_single,excluded=[3])
    return f_v(P,T,fo2,method=method)


######################################################
####################### A18 ###########################
def teos(T,P,V0,T0,dVdT,K0,Kp,KdP):
#    print("here")
#    print("Kp is",Kp)
#    print("K0 is",K0)
#    print("Kdp is",KdP)
    tmp = 1 + Kp + K0*KdP
    a = (1 + Kp)/tmp
    b = Kp/K0 - KdP/(1+Kp)
    c = tmp/(Kp**2 + Kp - K0*KdP)
    
    V0T = V0 + (T-T0)*dVdT
    V_V0T = 1-a*(1-(1+b*P)**(-c))  
    V = V0T*V_V0T
#    PV = V0T*P*(a*(b*P + 1)**(1-c)/(b*P*(1-c)) + (1-a))
    tmp2 = (1-a)*P + a/(b*(c-1)) - a*((1+b*P)**(1-c))/(b*(c-1)) # S10 in A19
    PV   = V0T*tmp2 # V is cm^3/mol, P is GPa, 
    return PV,V


## ox -0.10317460317460318, -Kp/K0 = -8/37
## re -0.21621621621621623, -Kp/K0 = -1.3/12.6
## What I fitted
# ox -0.10921177930813025
# re -0.1695887556539382
def fo2_A18(P,T,r,KdP_re =-0.21621621621621623, KdP_ox = -0.10317460317460318, method='earth'):  
    """
    
    Notes
    -------
    in other folders, this is wrong, np.log10(np.exp(1))*4*dPV/R/T
    also, from teos, to PV is wrong, 
    in this version, everything is corrects
    """
    Fe2 = (1-r)*planets.loc[method]['FeO']
    Fe3 = r*planets.loc[method]['FeO']    
    Si = planets.loc[method]['SiO2'];
    Mg = planets.loc[method]['MgO'];
    Al = planets.loc[method]['Al2O3'];
    Ca = planets.loc[method]['CaO'];
    K  = planets.loc[method]['K2O'];
    Na = planets.loc[method]['Na2O'];
    Ph = planets.loc[method]['P2O5'];    
    Ti = planets.loc[method]['TiO2'];
    data = A18_FeO_re
    PV_term_re,_ = teos(T,P,data['V0'],data['T0'],data['dVdT'],data['K0'],data['Kp'],KdP_re)
    data = A18_FeO_ox
    PV_term_ox,_ = teos(T,P,data['V0'],data['T0'],data['dVdT'],data['K0'],data['Kp'],KdP_ox)
    dPV= (PV_term_ox - PV_term_re)*1e3
#    print("dPV/R/T is",np.log10(np.exp(4*dPV/R/T)))
    x = r/(1-r)
#    other = 4*np.log10(x)-28144/T + 13.95 + 3905*Mg/T - 13359*Ca/T - 14858*Na/T - \
#         9805*K/T + 10906*Al/T + 110971*Ph/T - 11952*(Fe2 -Fe3)/T 
#    print("others are", other)
    log10fo2 = 4*np.log10(x)-28144/T + 13.95 + 3905*Mg/T - 13359*Ca/T - 14858*Na/T - \
         9805*K/T + 10906*Al/T + 110971*Ph/T - 11952*(Fe2 -Fe3)/T + \
         + np.log10(np.exp(1))*4*dPV/R/T
#    print("dPV is",np.log10(np.exp(1))*4*dPV/R/T)
#    print("PV/RT is",dPV/R/T)
#    print("W term is",3905*Mg/T - 13359*Ca/T - 14858*Na/T - \
#         9805*K/T + 10906*Al/T + 110971*Ph/T - 11952*(Fe2 -Fe3)/T)
#    print("dG term is",-28144/T + 13.95)
    return log10fo2

def X_A18_single(P,T,fo2,method='earth'):
    func = lambda r: fo2_A18(P,T,r,method=method) - fo2
    try:
 #       out = opt.brentq(func, 0.0005, 0.1) # for iw
        out = opt.brentq(func, 0.0005, .99999999) # for iw
        return out
    except:
        print("r =  0.0005, func is", fo2_A18(P,T,0.0005,method=method) - fo2)
        print("r =  .99, func is", fo2_A18(P,T,.5,method=method) - fo2)    

def X_A18(P,T,fo2,method='earth'):
    """
    calculate XFe3/Xtot as a function at fixed fo2, as a function of P
    -----
    input  list
    fo2 -> absolute fo2
    r   -> (Fe3/XFe)
    """

    f_v = np.vectorize(X_A18_single,excluded=[3])
    return f_v(P,T,fo2,method=method)


