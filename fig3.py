#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 00:11:29 2019
run fig3, then run fig Sx
@author: jiedeng
"""

# Earth: IW-2, 0-60 GPa
# Mars: IW-1.5, 0-15 GPa
# Moon: IW-2, 0-5 GPa
# Mercury: IW-3,0-5 GPa

## for Fe 25%
import src.tools2 as tl2
import src.tools as tl
import numpy as np
from scipy.interpolate import interp1d 
from src.geotherm_MgSiO3_melt import geotherm
import src.IW_buffer as iw
import matplotlib.pyplot as plt

def surface(P,Tgeo,dv=0,diw=-1.5,flag='Flase',method='Z17',PV_cal='old'):
    x_fe   = tl2.X_cal_vector(P,Tgeo,dv,np.log10(fiw)+diw,flag=False,method=method,PV_cal=PV_cal)
    size   = len(P)
    f_surf = np.zeros(size)
    for i in range(size):
        if i + 2 <  size:
            logf = tl2.fo2_cal_vector(P[:(i + 2)],Tgeo[:(i + 2)],dv[:(i + 2)],
                                      r = x_fe[i],flag=False,method=method,PV_cal=PV_cal)
        else:
            logf = tl2.fo2_cal_vector(P,Tgeo,dv,
                                      r = x_fe[i],flag=False,method=method,PV_cal=PV_cal)            
            
        f_surf[i] = logf[0]
    return f_surf,x_fe

def extrapolate_f(P,fs,p = [10,14]):
    f_fo2 = interp1d(P,fs)
    Pm = np.linspace(p[0],p[1],80)
    return Pm,f_fo2(Pm) - np.log10(fiw[0])

def extrapolate_r(P,r,p = [10,14]):
    f_r = interp1d(P,r)
    Pm = np.linspace(p[0],p[1],80)
    return Pm,f_r(Pm) 


Tm0          = 2100
P    = np.linspace(1e5,140e9,100)/1e9
Tgeo = geotherm(P*1e9,Tm0=Tm0)
_, _,_, _, dv_fe_25   = tl.cal_dV_this_study(P,Tgeo,name='Fe_25',flag=False)
_, _,_, _, dv_fe_12p5 = tl.cal_dV_this_study(P,Tgeo,name='Fe_12p5',flag=False)

#fiw = iw.f3_vector(P,Tgeo)
fiw = iw.f2_vector(P,Tgeo)
fiw_cold = fiw

fs_25_cold,r_25_cold                 = surface(P,Tgeo,dv=dv_fe_25,diw=-1.5,flag='False',method='mars',PV_cal = '25_cold')
fs_12p5_earth_cold,r_12p5_earth_cold = surface(P,Tgeo,dv=dv_fe_12p5,diw=-2,flag='False',method='earth',PV_cal = '12p5_cold')
fs_12p5_moon_cold, r_12p5_moon_cold  = surface(P,Tgeo,dv=dv_fe_12p5,diw=-2,flag='False',method='moon',PV_cal = '12p5_cold')

p_mars_cold,f_mars_cold   = extrapolate_f(P,fs_25_cold,[10,14]);
p_earth_cold,f_earth_cold = extrapolate_f(P,fs_12p5_earth_cold,[25,90]);
p_moon_cold,f_moon_cold   = extrapolate_f(P,fs_12p5_moon_cold,[1,7])

p_mars_cold,r_mars_cold   = extrapolate_r(P,r_25_cold,[10,14]);
p_earth_cold,r_earth_cold = extrapolate_r(P,r_12p5_earth_cold,[25,90]);
p_moon_cold,r_moon_cold   = extrapolate_r(P,r_12p5_moon_cold,[1,7])

Tm0          = 2500
#P    = np.linspace(1e5,80e9,100)/1e9
Tgeo = geotherm(P*1e9,Tm0=Tm0)
_, _,_, _, dv_fe_25 = tl.cal_dV_this_study(P,Tgeo,name='Fe_25',flag=False)
_, _,_, _, dv_fe_12p5 = tl.cal_dV_this_study(P,Tgeo,name='Fe_12p5',flag=False)

fiw = iw.f2_vector(P,Tgeo)
fiw_hot = fiw

fs_25_hot,r_25_hot                 = surface(P,Tgeo,dv=dv_fe_25,diw=-1.5,flag='False',method='mars',PV_cal = '25_hot')
fs_12p5_earth_hot,r_12p5_earth_hot = surface(P,Tgeo,dv=dv_fe_12p5,diw=-2,flag='False',method='earth',PV_cal = '12p5_hot')
fs_12p5_moon_hot, r_12p5_moon_hot  = surface(P,Tgeo,dv=dv_fe_12p5,diw=-2,flag='False',method='moon',PV_cal = '12p5_hot')

p_mars_hot,f_mars_hot   = extrapolate_f(P,fs_25_hot,[10,14]);
p_earth_hot,f_earth_hot = extrapolate_f(P,fs_12p5_earth_hot,[25,90]);
p_moon_hot,f_moon_hot   = extrapolate_f(P,fs_12p5_moon_hot,[1,7])

p_mars_hot,r_mars_hot   = extrapolate_r(P,r_25_hot,[10,14]);
p_earth_hot,r_earth_hot = extrapolate_r(P,r_12p5_earth_hot,[25,90]);
p_moon_hot,r_moon_hot   = extrapolate_r(P,r_12p5_moon_hot,[1,7])

alpha = 0.6

fig,ax = plt.subplots(1,2,figsize=(10,6),sharey=True)

from vatic.plots import plot_tool
plot_tool.load_default_setting()
fig.subplots_adjust(hspace=0.2,wspace=0.1)
#plot cold geotherm
#ax[0].plot(fs_25_cold - np.log10(fiw_cold[0]),P,'r-',alpha=alpha,label='')
ax[0].plot(f_mars_cold,p_mars_cold,'r-',linewidth=3)

#ax[0].plot(fs_12p5_earth_cold - np.log10(fiw_cold[0]),P,'b-',alpha=alpha,label='')
ax[0].plot(f_earth_cold,p_earth_cold,'b-',linewidth=3)

#ax[0].plot(fs_12p5_moon_cold - np.log10(fiw_cold[0]),P,'g-',alpha=alpha,label='')
ax[0].plot(f_moon_cold,p_moon_cold,'g-',linewidth=3)

#plot hot geotherm
#ax[0].plot(fs_25_hot - np.log10(fiw_hot[0]),P,'r--',alpha=alpha,label='')
#ax[0].plot(f_mars_hot,p_mars_hot,'r--',linewidth=3)
#
#ax[0].plot(fs_12p5_earth_hot - np.log10(fiw_hot[0]),P,'b--',alpha=alpha,label='')
#ax[0].plot(f_earth_hot,p_earth_hot,'b--',linewidth=3)
#
#ax[0].plot(fs_12p5_moon_hot - np.log10(fiw_hot[0]),P,'g--',alpha=alpha,label='')
#ax[0].plot(f_moon_hot,p_moon_hot,'g--',linewidth=3)

ax[0].grid(False)
ax[0].set_ylabel('MO base pressure (GPa)',fontsize=14)
ax[0].set_xlabel('Surface redox (' + r'$\mathrm{log }f_{\mathrm{O}_2}$'+  ' - IW)',fontsize=14)
ax[0].set_ylim([0,100]);
ax[0].set_xlim([-3,5])
xrange = list(range(-3,6))
ax[0].set_xticks(xrange)
ax[0].invert_yaxis()

ax[1].grid(False)
ax[1].set_xlabel(r'$\mathrm{Fe}^{3+} /\Sigma \mathrm{Fe}$',fontsize=14)
ax[1].set_ylim([0,100]);ax[1].set_xlim([0,0.042])
ax[1].invert_yaxis()

### plot Fe3/Fe
#ax[1].plot(r_25_cold,P,'r',label='Mars cold',alpha=alpha)
ax[1].plot(r_mars_cold,p_mars_cold,'r-',linewidth=3,label='')

#ax[1].plot(r_12p5_earth_cold,P,'b',label='Earth cold',alpha=alpha)
ax[1].plot(r_earth_cold,p_earth_cold,'b-',linewidth=3,label='')

#ax[1].plot(r_12p5_moon_cold,P,'g-',label='Moon cold',alpha=alpha)
ax[1].plot(r_moon_cold,p_moon_cold,'g-',linewidth=3,label='')

#ax[1].plot(r_25_hot,P,'r--',label='Mars hot',alpha=alpha)
#ax[1].plot(r_mars_hot,p_mars_hot,'r--',linewidth=3,label='')
#
#ax[1].plot(r_12p5_earth_hot,P,'b--',label='Earth hot',alpha=alpha)
#ax[1].plot(r_earth_hot,p_earth_hot,'b--',linewidth=3,label='')
#
#ax[1].plot(r_12p5_moon_hot,P,'g--',label='Moon hot',alpha=alpha)
#ax[1].plot(r_moon_hot,p_moon_hot,'g--',linewidth=3,label='')

ax[0].plot(np.arange(2.5,5.5),np.zeros_like(np.arange(2.5,5.5))-3.5,'b-',linewidth = 5,clip_on = False,label='Earth')
ax[0].plot(np.arange(-1,2),np.zeros_like(np.arange(-1,2))-2.5,'r-',linewidth = 5,clip_on = False,label='Mars')
ax[0].plot(np.arange(-2,1),np.zeros_like(np.arange(-2,1))-0.5,'g-',linewidth = 5,clip_on = False,label='Moon')

ax[1].plot(np.linspace(0.01,0.04,4),np.zeros_like(np.linspace(0.01,0.04,4))-3.5,'b-',linewidth = 5,clip_on = False,label='Earth')
#ax[1].plot(np.linspace(0.002,0.0004,4),np.zeros_like(np.linspace(0.002,0.0004,4))-2.5,'r-',linewidth = 5,clip_on = False,label='Mars')

#ax[1].plot(np.arange(-1,2),np.zeros_like(np.arange(-1,2))-1.5,'r-',linewidth = 5,clip_on = False,label='Mars')
#ax[1].plot(np.arange(-2,1),np.zeros_like(np.arange(-2,1))-0.5,'g-',linewidth = 5,clip_on = False,label='Moon')
#ax[0].plot([0,1],[-1,-2],'k-',linewidth = 3,label='cold')
#ax[0].plot([0,1],[-1,-2],'k--',alpha=0.7,label='hot')
ax[0].legend(fontsize=14)

ax[0].text(-2.8,97,'a',fontsize=16, fontweight='bold')
ax[1].text(0.0016,97,'b',fontsize=16, fontweight='bold')
plot_tool.show_minor_ticks(ax)
plot_tool.set_major_axis_font(ax[0],14)
plot_tool.set_major_axis_font(ax[1],14)


ax_extra_y = ax[1].twinx()
ax_extra_y.set_ylabel('MO base depth (km)',fontsize=14)
yticks = np.linspace(0,100,6)

import burnman
prem = burnman.seismic.PREM()
y_depth = prem.depth(yticks*1e9)   ### need revise the burnman code from xx to xx.any()
y_depth = np.round(y_depth/1e3)
ax_extra_y.set_yticks(yticks)
ax_extra_y.set_yticklabels([str(int(i)) for i in y_depth],fontsize=14)
ax_extra_y.invert_yaxis()
ax_extra_y.minorticks_on()


if False:
    fig.savefig('out/fig3_only_cold_2.pdf',bbox_inches='tight')
