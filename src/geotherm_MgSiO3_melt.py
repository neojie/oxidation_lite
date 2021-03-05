#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 00:18:02 2019

@author: jiedeng
"""
import numpy as np
from scipy.interpolate import interp1d
from burnman import Mineral

        
class mgsio3(Mineral):
    """
    eos of mgsio3 melt
    ref: Liebske and Frost, 2012, EPSL
    """
    def __init__(self):
        self.params = {
            'equation_of_state': 'lf12l',
            'T_0': 1773,
            'V_0': 37.2e-6,  
            'K_0': 27.3e9,
            'Kprime_0': 5.7,
            'molar_mass': .1,
            'n': 5, 
            'grueneisen_0': 0.6,
            'grueneisen_prime': -1.24,  # dG/dV
            'Cv': 173,
            'S_0': 333,
            'F_0':-1710300}

        Mineral.__init__(self)

def geotherm(P,Tm0=2500):
#    pv_l = minerals.LF_2012.mgsio3()
    pv_l = mgsio3()
    V0  = pv_l.evaluate(['molar_volume'],[0],[Tm0])[0];
    V   = np.linspace(V0/3,V0,100);
    tmp = np.log(V/V0)*(pv_l.params['grueneisen_prime'] - pv_l.params['grueneisen_0']) - \
          pv_l.params['grueneisen_prime']/V0*(V-V0)
    Tm  = Tm0*np.exp(tmp);
    Pm  = np.zeros(V.shape);
    for i in range(len(V)):
        Pm[i] = pv_l.method.pressure(Tm[i], V[i], pv_l.params)   

    f   = interp1d(Pm[:,0],Tm[:,0])

    Tin = f(P)
    return Tin






