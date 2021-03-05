#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 22:43:56 2018

@author: jiedeng
"""
from scipy.interpolate import interp1d
import numpy as np
import scipy.constants as con

def linfit(x,k,y0):
    """
    The polyfit does not give std 
    """
    y = x*k+y0
    return y

def P2V(P0,V0,Pint):
    """
    interpolate pressure to volume
    """
    f_int = interp1d(P0,V0)
    V_int = f_int(Pint)
    return V_int

def eval_alpha(V,T,P):
    """    
    V:volume
    T:temperature
    P:pressure  
    assumtions: 1 - V vs. T at const P is a straight line
                2 - average V when dV/V    
    """
    # Nint control the interpolate size
    Nint = 100
    row,col = np.shape(P)
    
    if (row is np.shape(T)[0]) and (col is np.shape(V)[0]):
        print("input data format correct")
    else:
        print("Wrong, data dimension not match")
        
    if V[0]>V[-1]:
        Pmin = np.max(P[:,0])
        Pmax = np.min(P[:,-1])
    else:
        Pmin = np.max(P[:,-1])
        Pmax = np.min(P[:,0])    
        
    Pint = np.linspace(Pmin,Pmax,Nint)   
    #print(Pint)
    Vint  = np.zeros([row,Nint])
    alpha = np.zeros([Nint])
    
    for i in np.arange(row):
        Vint[i] = P2V(P[i],V,Pint) 
               
    for i in np.arange(Nint):
        dVdT  = np.polyfit(T,Vint[:,i],1)[0]
#        print(dVdT)
#        print(Vint[:,i])
        alpha[i] = dVdT/np.mean(Vint[:,i])
    return alpha,Pint,Vint

def eval_KT(P,V,M= 55.85+16*2+1,P_tar=0):
    """
    sf: pressure
    sx: volume
    parameters: array of eos para
    P: target pressure
    M:Molar mass of FeO2H g/mol
    """
    # 
    # direct evaluation of KT
    #Avogadro number
    Na  = con.N_A   
    # python interpolation -> volume per formula unit at target pressure P GPa A^3
    # Density kg/m^3
    rho      = M/(V*1e-30*Na)/1e3
    dPdV     = np.gradient(P, V)   #dP/dV
    dPdV_int = interp1d(P,dPdV)
    KT = -dPdV_int(P)*V
    return rho, KT

def bm(vol,T,K0,Kp,Kdp,V0,P0,T0,a,b,c):
    """
    1st version
    """
    B1 = K0*Kdp+(Kp-4.)*(Kp-5.)+59./9. #these . are critical
    B2 = 3*K0*Kdp+(Kp-4)*(3*Kp-13)+129./9.
    B3 = 3*K0*Kdp+(Kp-4)*(3*Kp-11)+105./9.
    B4 = K0*Kdp+(Kp-4)*(Kp-3)+35./9.   
    f  = vol/V0  
    BH = (a - b*f +c*(f**2))/1000.
    P = 9./16.*K0*(-B1*(f**(-5./3.))+B2*(f**(-7./3.))\
            -B3*(f**(-3.))+B4*(f**(-11./3.)))+P0+BH*(T-T0)
    return P

def bm4(vol,K0,Kp,Kdp,V0,P0,T=0):
    """
    from bm, but remove the thermal part
    """
    B1 = K0*Kdp+(Kp-4.)*(Kp-5.)+59./9. #these . are critical
    B2 = 3*K0*Kdp+(Kp-4)*(3*Kp-13)+129./9.
    B3 = 3*K0*Kdp+(Kp-4)*(3*Kp-11)+105./9.
    B4 = K0*Kdp+(Kp-4)*(Kp-3)+35./9.   
    f  = vol/V0  
    Pc = 9./16.*K0*(-B1*(f**(-5./3.))+B2*(f**(-7./3.))\
            -B3*(f**(-3.))+B4*(f**(-11./3.)))+P0
    return Pc

def BM3(vol,K0,Kp,V0,P0):
    V = vol
    a1 = 3*V0*P0
    a2 = 3*V0*(3*K0-5*P0)/2
    a3 = V0*(9*K0*Kp - 36*K0 + 35*P0)/2
    f = 1/2*((V0/V)**(2/3)-1) 
    dfdV = - (V0**(2/3))/3/(V**(5/3)) 
    P_jd3 = (a1 + 2*a2*f + 3*a3*(f**2))*(-dfdV)
    return P_jd3

def BM4_single_v(vol,K0,Kp,Kdp,V0,P0):
    """
    Note when Kdp ==0, BM4 does not go degenerate to BM3, because Kdp and Kp, K0
    are strongly correlated with each other!!
    """
    V = vol
    if vol >0: #if vol <= V0: this one does not help
        a1 = 3*V0*P0
        a2 = 3*V0*(3*K0-5*P0)/2
        a3 = V0*(9*K0*Kp - 36*K0 + 35*P0)/2
        
        a4 = 3*V0*(9*(K0**2)*Kdp + 9*K0*(Kp**2) 
                   - 63*K0*Kp + 143*K0 - 105*P0)/8.0
        f = 1/2*((V0/V)**(2/3)-1) 
        dfdV = - (V0**(2/3))/3/(V**(5/3)) 
        P_jd4 = (a1 + 2*a2*f + 3*a3*(f**2) + 4*a4*(f**3))*(-dfdV)
    else:
        P_jd4 = 0
    return P_jd4

def BM4(vol,K0,Kp,Kdp,V0,P0):
    """
    Oct10,2018, I add a4 1/8!!!, did not notice before!
    """
    f_v = np.vectorize(BM4_single_v, excluded=list(range(1,6)))
    return f_v(vol,K0,Kp,Kdp,V0,P0)

def TH_single_v(vol,V0,T,T0,a,b,c):
    f  = vol/V0 
    BH = (a - b*f +c*(f**2))/1000.
    Pth = BH*(T-T0)
    return Pth

def TH(vol,V0,T,T0,a,b,c):
    th_v = np.vectorize(TH_single_v, excluded=[1,3,4,5,6])
    return th_v(vol,V0,T,T0,a,b,c)

def TH4(vol,V0,T,T0,a,b,c,d):
    f  = vol/V0 
    BH = (a - b*f +c*(f**2)+d*(f**3))/1000.
    Pth = BH*(T-T0)
    return Pth

def TH_LF12(vol,V0,T,T0,a,b,c,d):
    f  = vol/V0 
    gamma = a + b*(f-1) + c*((f-1)**2)
    Cv    = d
    Pth = gamma/vol*Cv*(T-T0)
    return Pth
 
def BM3_TH(vol,K0,Kp,V0,P0,T,T0,a,b,c):
    Pc  = BM3(vol,K0,Kp,V0,P0)
    Pth = TH(vol,V0,T,T0,a,b,c)
    return Pc+Pth

def BM4_TH(vol,K0,Kp,Kdp,V0,P0,T,T0,a,b,c):
    Pc  = BM4(vol,K0,Kp,Kdp,V0,P0)
    Pth = TH(vol,V0,T,T0,a,b,c)
    return Pc+Pth

def BM4_TH4(vol,K0,Kp,Kdp,V0,P0,T,T0,a,b,c,d):
    Pc  = BM4(vol,K0,Kp,Kdp,V0,P0)
    Pth = TH4(vol,V0,T,T0,a,b,c,d)
    return Pc+Pth

def BM4_TH_LF12(vol,K0,Kp,Kdp,V0,P0,T,T0,a,b,c,d):
    Pc  = BM4(vol,K0,Kp,Kdp,V0,P0)
    Pth = TH_LF12(vol,V0,T,T0,a,b,c,d)
    return Pc+Pth

def BM3_alpha(vol,K0,Kp,V0,P0,T,T0,dKdT,a,b):
    """
    ref : Liu 2014
    """
    V0T = V0*np.exp((T-T0)*a + (T**2 - T0**2)/2*b)
#    print("exp is",(T-T0)*a + (T**2 - T0**2)/2*b)
#    print("V0T is",V0T)
    K0T = K0 + dKdT*(T-T0)
    PT = BM3(vol,K0T,Kp,V0T,P0)
    return PT

def BM4_alpha(vol,K0,Kp,Kdp,V0,P0,T,T0,dKdT,a,b):
    """
    ref : 
    dKdT assumes to be a constant, which is not the case for liquid
    for liquid
    dKdT = dKdT0 + (P-P0)*dKdT_slope      
    I may create a BM4_alpha_l, remove one freedom for liquid 
    and add one for dKdT
    """
    V0T = V0*np.exp((T-T0)*a + (T**2 - T0**2)/2*b)
    print("exp is",(T-T0)*a + (T**2 - T0**2)/2*b)
    print("V0T is",V0T)
    K0T = K0 + dKdT*(T-T0)
    print("K0T is",K0T)
    PT = BM4(vol,K0T,Kp,Kdp,V0T,P0)
    return PT


def Vinet(vol,a,b,c,d,e,f,V0,K0,Kp):
    """
    """