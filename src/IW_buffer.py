#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 20:56:38 2019

@author: jiedeng
"""

import numpy as np
#### f1,f2,f3 are IW buffer

def f1(P,T):
    """
    calculate oxygen fugacity of IW buffer, Fe0.94O-Fe solid

    Parameters
    ----------
    input:   T in K; P in GPa
    output:  f-> absolute oxygen fugacity
    
    
    Ref.  eqn S1 in Zhang et al., 2017; Huebner, 1971
    """
    f = -28777.89/T + 14.0572 - 2.039*np.log10(T) + 550*(P-0.0001)/T
    return 10**f

def f2(P,T):
    """
    calculate oxygen fugacity of IW buffer, FeO solid - Fe solid
    Ref. Campbell et al., 2009, coefficents taken from SI
    """
    a0 = 6.54106; a1 = 0.0012324
    b0 = -28163.6; b1 = 546.32; b2 = -1.13412; b3 = 0.0019274
    f = (a0+a1*P) + (b0+b1*P+b2*(P**2)+b3*(P**3))/T
    return 10**f

def f2_vector(P,T):
    """
    calculate oxygen fugacity of IW buffer, FeO solid - Fe solid
    Ref. Campbell et al., 2009, coefficents taken from SI
    """
    f2_v = np.vectorize(f2)
    return f2_v(P,T)


def f3(P,T):
    """
    calculate oxygen fugacity of IW buffer, Fe0.94O

    Parameters
    ----------
    input:   T in K; P in GPa
    output:  f-> absolute oxygen fugacity
    
    ref. Zhang et al., 2017
    """
    f_campbell  = f2(P,T)
    f_composite = np.log10(f1(P,T))*(1-0.2*P) + np.log10(f2(P,T))*0.2*P 
    return np.maximum(f_campbell,10**f_composite)

def f3_vector(P,T):
    f_v = np.vectorize(f3)
    return f_v(P,T)
    
def f4(P,T):
    """
    ref. ONeill (1988)
    """
    

def f5(P,T):
    """
    calculate oxygen fugacity of FeO liquid-Fe solid   
    ref. 
    """

def f6(P,T):
    """
    calculate oxygen fugacity of FeO in silicate melt- Fe solid   
    """    
def Ru_RuO2_single(p,t):
    tmp = -16953/t+17.98-2.66*np.log10(t)+526*p/t
    return 10**tmp

def Ru_RuO2(p,t):
    f_v2 = np.vectorize(Ru_RuO2_single)
    return f_v2(p,t)


def Ru_RuO2_2_single(P,T):
    """
    page 77 of Armstrong thesis

    """
    a0 = 7.783; a1 = -0.0113; a2 = 0.00209; a3 = -4.1e-5; 
    b0 = -13764; b1 = 593.4;  b2 = -4.08; 
    c0 = -1.05e6; c1 = -4511
    f = (a0+a1*P+a2*(P**2)+a3*(P**3)) + (b0+b1*P+b2*(P**2))/T + (c0+c1*P)/(T**2)
    return 10**f

def Ru_RuO2_2(P,T):
    """
    page 77 of Armstrong thesis

    """
    f_v3 = np.vectorize(Ru_RuO2_2_single)
    return f_v3(P,T)