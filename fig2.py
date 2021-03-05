#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:39:34 2019

@author: jiedeng
"""

import src.IW_buffer as iw
from src.geotherm_MgSiO3_melt import geotherm
import src.tools as tl
import src.tools2 as tl2
import numpy as np
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt

def extrapolate(P,material,pivot = 25):
    f_fo2 = interp1d(P,material)
    f_iw  = interp1d(P,np.log10(fiw))
    Pm = np.linspace(1e-4,pivot,80)
    return Pm,f_fo2(Pm) - f_iw(Pm),f_fo2(Pm)

R = 8.314
A3_to_cm3   = 1e-30*1e6*6.022e23
cm3_to_A3   = 1e30*1e-6/6.022e23

Pearth = 55
# should be 55

Tm0  = 2100
P    = np.linspace(1e5,60e9,100)/1e9  ## should have more than 50 bins
Tgeo = geotherm(P*1e9,Tm0=Tm0) 
#Tgeo = geotherm(55*1e9,Tm0=Tm0)

_, _,_, _, dv_fe_25   = tl.cal_dV_this_study(P,Tgeo,name='Fe_25',flag=False)
_, _,_, _, dv_fe_12p5 = tl.cal_dV_this_study(P,Tgeo,name='Fe_12p5',flag=False)

#fiw = iw.f3_vector(P,Tgeo)
fiw = iw.f2_vector(P,Tgeo)
fiw_cold = fiw

# Earth: IW-2, 25 GPa
# Mars: IW-1.5, 15 GPa
# Moon: IW-2, 5 GPa
# Mercury: IW-3,5 GPa

## for Fe 25%
x_fe_25_cold = tl2.X_cal_vector(P,Tgeo,dv_fe_25,np.log10(fiw)-1.5,flag=False,method='mars', PV_cal='25_cold')
f_25_cold    = interp1d(P,x_fe_25_cold);
r_25_cold    = f_25_cold(14)
logf_25_cold = tl2.fo2_cal_vector(P,Tgeo,dv_fe_25,r = r_25_cold,flag=False,method='mars', PV_cal='25_cold')

# for Fe 12.5%
x_12p5_earth_cold = tl2.X_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,flag=False,method='earth', PV_cal='12p5_cold')
x_12p5_moon_cold  = tl2.X_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,flag=False,method='moon',  PV_cal='12p5_cold')

f_12p5_earth_cold = interp1d(P,x_12p5_earth_cold);r_12p5_earth_cold = f_12p5_earth_cold(Pearth);
f_12p5_moon_cold  = interp1d(P,x_12p5_moon_cold);r_12p5_moon_cold = f_12p5_moon_cold(5);

logf_12p5_earth_cold = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_12p5_earth_cold,flag=False,method='earth', PV_cal='12p5_cold')
logf_12p5_moon_cold  = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_12p5_moon_cold,flag=False,method='moon', PV_cal='12p5_cold')

p_mars_cold,f_mars_cold,fo_mars_cold   = extrapolate(P,logf_25_cold,14);
p_earth_cold,f_earth_cold,fo_earth_cold = extrapolate(P,logf_12p5_earth_cold,Pearth);
p_moon_cold,f_moon_cold,fo_moon_cold   = extrapolate(P,logf_12p5_moon_cold,5);

### Z17
x_Z17_fe_25_cold = tl2.X_Z17(P,Tgeo,np.log10(fiw)-1.5)
f_Z17_25_cold    = interp1d(P,x_Z17_fe_25_cold);
r_Z17_25_cold    = f_Z17_25_cold(14)
logf_Z17_25_cold = tl2.fo2_Z17(P,Tgeo,r = r_Z17_25_cold)

x_Z17_12p5_cold = tl2.X_Z17(P,Tgeo,np.log10(fiw)-2)

f_Z17_12p5_cold = interp1d(P,x_Z17_12p5_cold);
r_Z17_12p5_earth_cold = f_Z17_12p5_cold(25);
r_Z17_12p5_moon_cold = f_Z17_12p5_cold(5);

logf_Z17_12p5_earth_cold = tl2.fo2_Z17(P,Tgeo,r = r_Z17_12p5_earth_cold)
logf_Z17_12p5_moon_cold = tl2.fo2_Z17(P,Tgeo,r = r_Z17_12p5_moon_cold)

p_Z17_mars_cold,f_Z17_mars_cold,fo_Z17_mars_cold   = extrapolate(P,logf_Z17_25_cold,14);
p_Z17_earth_cold,f_Z17_earth_cold,fo_Z17_earth_cold = extrapolate(P,logf_Z17_12p5_earth_cold,25);
p_Z17_moon_cold,f_Z17_moon_cold,fo_Z17_moon_cold   = extrapolate(P,logf_Z17_12p5_moon_cold,5);

### O16
x_O06_fe_25_cold = tl2.X_O06(P,Tgeo,np.log10(fiw)-1.5,method='mars')
f_O06_25_cold    = interp1d(P,x_O06_fe_25_cold);
r_O06_25_cold    = f_O06_25_cold(14)
logf_O06_25_cold = tl2.fo2_O06(P,Tgeo,r = r_O06_25_cold,method='mars')

x_O06_12p5_earth_cold = tl2.X_O06(P,Tgeo,np.log10(fiw)-2,method='earth')
x_O06_12p5_moon_cold = tl2.X_O06(P,Tgeo,np.log10(fiw)-2,method='moon')

f_O06_12p5_earth_cold = interp1d(P,x_O06_12p5_earth_cold);
f_O06_12p5_moon_cold  = interp1d(P,x_O06_12p5_moon_cold);

r_O06_12p5_earth_cold = f_O06_12p5_earth_cold(Pearth);
r_O06_12p5_moon_cold  = f_O06_12p5_moon_cold(5);

logf_O06_12p5_earth_cold = tl2.fo2_O06(P,Tgeo,r = r_O06_12p5_earth_cold,method='earth')
logf_O06_12p5_moon_cold = tl2.fo2_O06(P,Tgeo,r = r_O06_12p5_moon_cold,method='moon')

p_O06_mars_cold,f_O06_mars_cold,fo_O06_mars_cold   = extrapolate(P,logf_O06_25_cold,14);
p_O06_earth_cold,f_O06_earth_cold,fo_O06_earth_cold = extrapolate(P,logf_O06_12p5_earth_cold,Pearth);
p_O06_moon_cold,f_O06_moon_cold,fo_O06_moon_cold   = extrapolate(P,logf_O06_12p5_moon_cold,5);

### A18
x_A18_fe_25_cold = tl2.X_A18(P,Tgeo,np.log10(fiw)-1.5,method='mars')
f_A18_25_cold    = interp1d(P,x_A18_fe_25_cold);
r_A18_25_cold    = f_A18_25_cold(14)
logf_A18_25_cold = tl2.fo2_A18(P,Tgeo,r = r_A18_25_cold,method='mars')

x_A18_12p5_earth_cold = tl2.X_A18(P,Tgeo,np.log10(fiw)-2,method='earth')
x_A18_12p5_moon_cold  = tl2.X_A18(P,Tgeo,np.log10(fiw)-2,method='moon')

f_A18_12p5_earth_cold = interp1d(P,x_A18_12p5_earth_cold);
f_A18_12p5_moon_cold  = interp1d(P,x_A18_12p5_moon_cold);

r_A18_12p5_earth_cold = f_A18_12p5_earth_cold(25);
r_A18_12p5_moon_cold  = f_A18_12p5_moon_cold(5);

logf_A18_12p5_earth_cold = tl2.fo2_A18(P,Tgeo,r = r_A18_12p5_earth_cold,method='earth')
logf_A18_12p5_moon_cold = tl2.fo2_A18(P,Tgeo,r = r_A18_12p5_moon_cold,method='moon')

p_A18_mars_cold,f_A18_mars_cold,fo_A18_mars_cold   = extrapolate(P,logf_A18_25_cold,14);
p_A18_earth_cold,f_A18_earth_cold,fo_A18_earth_cold = extrapolate(P,logf_A18_12p5_earth_cold,25);
p_A18_moon_cold,f_A18_moon_cold,fo_A18_moon_cold   = extrapolate(P,logf_A18_12p5_moon_cold,5);



fig,ax = plt.subplots(1,1,figsize=(3.8,4.5),sharey=True)

from vatic.plots import plot_tool
plot_tool.load_default_setting()
ax.plot(f_mars_cold,p_mars_cold,'r-',linewidth = 3,label='')
ax.plot(f_earth_cold,p_earth_cold,'b-',linewidth = 3,label='')
ax.plot(f_moon_cold,p_moon_cold,'g-',linewidth = 3,label='')

ax.plot(f_Z17_mars_cold,p_Z17_mars_cold,'r--',alpha=0.7,label='')
ax.plot(f_Z17_earth_cold,p_Z17_earth_cold,'b--',alpha=0.7,label='')
ax.plot(f_Z17_moon_cold,p_Z17_moon_cold,'g--',alpha=0.7,label='')

#ax.plot(f_O06_mars_cold,p_O06_mars_cold,'r:',alpha=0.7,label='')
#ax.plot(f_O06_earth_cold,p_O06_earth_cold,'b:',alpha=0.7,label='')
ax.plot(f_O06_moon_cold,p_O06_moon_cold,'g:',alpha=0.7,label='')

ax.plot(f_A18_mars_cold,p_A18_mars_cold,'r-.',alpha=0.7,label='')
ax.plot(f_A18_earth_cold,p_A18_earth_cold,'b-.',alpha=0.7,label='')
ax.plot(f_A18_moon_cold,p_A18_moon_cold,'g-.',alpha=0.7,label='')

ax.plot(np.arange(2.5,5.5),np.zeros_like(np.arange(2.5,5.5))-2.5,'b-',linewidth = 5,clip_on = False,label='Earth')
ax.plot(np.arange(-1,2),np.zeros_like(np.arange(-1,2))-1.5,'r-',linewidth = 5,clip_on = False,label='Mars')
ax.plot(np.arange(-2,1),np.zeros_like(np.arange(-2,1))-0.5,'g-',linewidth = 5,clip_on = False,label='Moon')

ax.plot([0,1],[-1,-2],'k-',linewidth = 3,label='This study')
ax.plot([0,1],[-1,-2],'k-.',alpha=0.7,label='Armstrong 2018')
ax.plot([0,1],[-1,-2],'k--',alpha=0.7,label='Zhang et al., 2017')
ax.plot([0,1],[-1,-2],'k:',alpha=0.7,label='O\'Neill et al,, 2006')

#ax.set_xlabel('Redox within MO (' + r'$\mathrm{log } f_{O_{2}}$'+  ' - IW)',fontsize=14)
ax.set_xlabel('Redox within MO (' + r'$\mathrm{log }$'+ '$f_{O_{2}}$'+  ' - IW)',fontsize=14,fontname = 'Times')

ax.legend()
ax.set_ylim([0,Pearth])
ax.invert_yaxis()

ax.set_ylabel('Pressure within MO (GPa)',fontsize=14)
ax.set_ylim([0,Pearth])
ax.set_xlim([-5,5])
ax.invert_yaxis()
yticks = np.linspace(0,50,6)
yticks = np.append(yticks,[55])
ax.set_yticks(yticks)
ax.set_yticklabels([str(int(i)) for i in yticks])

xticks = np.arange(-5,6,2)
ax.set_xticks(xticks)
ax.set_xticklabels([str(int(i)) for i in xticks])



ax_extra_y = ax.twinx()
ax_extra_y.set_ylabel('Depth within MO (km)',fontsize=14)

import burnman
prem = burnman.seismic.PREM()
y_depth = prem.depth(yticks*1e9)   ### need revise the burnman code from xx to xx.any()
y_depth = np.round(y_depth/1e3)
ax_extra_y.set_yticks(yticks)
ax_extra_y.set_yticklabels([str(int(i)) for i in y_depth],fontsize=13.5)
ax_extra_y.invert_yaxis()
ax_extra_y.minorticks_on()


plot_tool.show_minor_ticks(ax)
plot_tool.set_major_axis_font(ax,13.5)

if False:
    fig.savefig('out/fig2_.pdf',bbox_inches='tight')




