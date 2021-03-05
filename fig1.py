#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 11:40:17 2019
generate fig 1
@author: jiedeng
"""

import src.tools as tl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plot_tool


P = np.linspace(0,140e9,141)/1e9

T = np.ones(P.shape)*2000
_, _,v_fe_25_re_2k, v_fe_25_ox_2k, dv_fe_25_2k = tl.cal_dV_this_study(P,T,name='Fe_25',flag=False)
_, _,v_fe_12p5_re_2k, v_fe_12p5_ox_2k, dv_fe_12p5_2k = tl.cal_dV_this_study(P,T,name='Fe_12p5',flag=False)

#T = np.ones(P.shape)*2500
#_, _,v_fe_25_re_25, v_fe_25_ox_25, dv_fe_25_25 = tl.cal_dV_this_study(P,T,name='Fe_25',flag=False)
#_, _,v_fe_12p5_re_25, v_fe_12p5_ox_25, dv_fe_12p5_25 = tl.cal_dV_this_study(P,T,name='Fe_12p5',flag=False)

T = np.ones(P.shape)*3000
_, _,v_feo_re_3k, v_feo_ox_3k, dv_feo_3k         = tl.cal_dV_this_study(P,T,name='FeO',flag=False)
_, _,v_fe_25_re_3k, v_fe_25_ox_3k, dv_fe_25_3k   = tl.cal_dV_this_study(P,T,name='Fe_25',flag=False)
_, _,v_fe_12p5_re_3k, v_fe_12p5_ox_3k, dv_fe_12p5_3k = tl.cal_dV_this_study(P,T,name='Fe_12p5',flag=False)

T = np.ones(P.shape)*4000
_, _,v_feo_re_4k, v_feo_ox_4k, dv_feo_4k         = tl.cal_dV_this_study(P,T,name='FeO',flag=False)
_, _,v_fe_25_re_4k, v_fe_25_ox_4k, dv_fe_25_4k   = tl.cal_dV_this_study(P,T,name='Fe_25',flag=False)
_, _,v_fe_12p5_re_4k, v_fe_12p5_ox_4k, dv_fe_12p5_4k = tl.cal_dV_this_study(P,T,name='Fe_12p5',flag=False)
#_, _,v_fe_6p25_re_4k, v_fe_6p25_ox_4k, dv_fe_6p25_4k = tl.cal_dV_this_study(P,T,name='Fe_6p25',flag=False)

xlsx = 'db/eos.xlsx'

def load_raw_data(sheet_name = 'melt6_FeO',skiprow = 1,nrow=18):                 
    FeO = pd.read_excel(xlsx,sheet_name=sheet_name,
                          usecols = list(range(11)), 
                          index_col = 0,
                          skiprows = list(range(skiprow)), 
                          nrows = nrow,).dropna()
    return FeO

fe_25_re = load_raw_data(sheet_name = 'melt_25',skiprow = 1,nrow=23)
fe_25_ox = load_raw_data(sheet_name = 'melt_25',skiprow = 23,nrow=26)
fe_12p5_re = load_raw_data(sheet_name = 'melt3_12p5',skiprow = 1,nrow=25)
fe_12p5_ox = load_raw_data(sheet_name = 'melt3_12p5',skiprow = 27,nrow=28)


Psel = 'P(GPa)'

####### plot  ######
fig,ax = plt.subplots(1,2,figsize=(11,4.5),sharex=True)
#plot_tool.load_default_setting()

fig.subplots_adjust(hspace=0.2,wspace=0.1)
left, bottom, width, height = [0.233, 0.4, 0.24, 0.4]
inset0 = fig.add_axes([left, bottom, width, height])
left, bottom, width, height = [0.641, 0.4, 0.24, 0.4]
inset1 = fig.add_axes([left, bottom, width, height])

ax[0].plot(P[:21],dv_fe_12p5_2k[:21],'g-',label= '2000 K')
ax[0].plot(P[:101],dv_fe_12p5_3k[:101],'b-',label= '3000 K')
ax[0].plot(P,dv_fe_12p5_4k,'r-',label= '4000 K')

re_name = '16 ' +r'$\mathrm{Mg}_{0.875}\mathrm{Fe}_{0.125}^{2+}\mathrm{Si}\mathrm{O}_{3}$';
ox_name = '16 ' +r'$\mathrm{Mg}_{0.875}\mathrm{Fe}_{0.125}^{3+}\mathrm{Si}\mathrm{O}_{3.0625}$';

inset0.plot(0,-100,'ko',markersize=4,label= re_name)
inset0.plot(0,-100,'ko',markerfacecolor='none',label= ox_name)


inset0.plot(P[:21],v_fe_12p5_re_2k[:21],'g-',label= ' 2000 K')
inset0.plot(P[:21],v_fe_12p5_ox_2k[:21],'g--')
inset0.plot(P[:101],v_fe_12p5_re_3k[:101],'b-',label= ' 3000 K')
inset0.plot(P[:101],v_fe_12p5_ox_3k[:101],'b--')
inset0.plot(P,v_fe_12p5_re_4k,'r-',label= ' 4000 K')
inset0.plot(P,v_fe_12p5_ox_4k,'r--')

inset0.plot(fe_12p5_re.loc[2000][Psel],fe_12p5_re.loc[2000]['V(A3)'],'go',markersize=4)
inset0.plot(fe_12p5_ox.loc[2000][Psel],fe_12p5_ox.loc[2000]['V(A3)'],'go',markerfacecolor='none',label= '')
inset0.plot(fe_12p5_re.loc[3000][Psel],fe_12p5_re.loc[3000]['V(A3)'],'bo',markersize=4)
inset0.plot(fe_12p5_ox.loc[3000][Psel],fe_12p5_ox.loc[3000]['V(A3)'],'bo',markerfacecolor='none',label= '')
inset0.plot(fe_12p5_re.loc[4000][Psel],fe_12p5_re.loc[4000]['V(A3)'],'ro',markersize=4)
inset0.plot(fe_12p5_ox.loc[4000][Psel],fe_12p5_ox.loc[4000]['V(A3)'],'ro',markerfacecolor='none',label= '')

inset0.set_xlim([0,140]);inset0.set_ylim([500,1400])
inset0_y = [ 500.,  725.,  950., 1175., 1400.]
inset0.set_yticks(inset0_y)
inset0.set_xlabel('P (GPa)', fontsize=13, fontname="Arial")
inset0.legend(frameon=False,fontsize=10,loc='best')

xtick = np.linspace(0,140,8)

inset1.plot(P[:21],v_fe_25_re_2k[:21],'g-')
inset1.plot(P[:101],v_fe_25_re_3k[:101],'b-')
inset1.plot(P,v_fe_25_re_4k,'r-')

inset1.plot(P[:21],v_fe_25_ox_2k[:21],'g--')
inset1.plot(P[:101],v_fe_25_ox_3k[:101],'b--')
inset1.plot(P,v_fe_25_ox_4k,'r--')

re_name = '16 ' + r'$\mathrm{Mg}_{0.75}\mathrm{Fe}_{0.25}^{2+}\mathrm{Si}\mathrm{O}_{3}$';
ox_name = '16 ' + r'$\mathrm{Mg}_{0.75}\mathrm{Fe}_{0.25}^{3+}\mathrm{Si}\mathrm{O}_{3.125}$';
inset1.plot(0,-100,'ko',markersize=4,label= re_name)
inset1.plot(0,-100,'ko',markerfacecolor='none',label= ox_name)

inset1.plot(fe_25_re.loc[2000][Psel],fe_25_re.loc[2000]['V(A3)'],'go',markersize=4,label= '')
inset1.plot(fe_25_re.loc[3000][Psel],fe_25_re.loc[3000]['V(A3)'],'bo',markersize=4,label= '')
inset1.plot(fe_25_re.loc[4000][Psel],fe_25_re.loc[4000]['V(A3)'],'ro',markersize=4,label= '')

inset1.plot(fe_25_ox.loc[2000][Psel],fe_25_ox.loc[2000]['V(A3)'],'go',markerfacecolor='none',label= '')
inset1.plot(fe_25_ox.loc[3000][Psel],fe_25_ox.loc[3000]['V(A3)'],'bo',markerfacecolor='none',label= '')
inset1.plot(fe_25_ox.loc[4000][Psel],fe_25_ox.loc[4000]['V(A3)'],'ro',markerfacecolor='none',label= '')
inset1.legend(frameon=False,fontsize=10)

inset1.set_xlim([0,140]);
inset1.set_ylim([500,1400]);
inset1_y = [ 500.,  725.,  950., 1175., 1400.]
inset1.set_yticks(inset1_y)
inset1.set_xlabel('P (GPa)', fontsize=13, fontname="Arial")

ax[1].plot(P[:21],dv_fe_25_2k[:21],'g-')
ax[1].plot(P[:101],dv_fe_25_3k[:101],'b-')
ax[1].plot(P,dv_fe_25_4k,'r-')

inset0.set_xticks(xtick)
inset1.set_xticks(xtick)

plot_tool.set_major_axis_font(inset0,12.5)
plot_tool.set_major_axis_font(inset1,12.5)

inset0.set_ylabel("Volume"+" ("+r"$\AA^3$"+")",fontsize=13)
inset1.set_ylabel("Volume"+" ("+r"$\AA^3$"+")",fontsize=13)

ax[0].set_ylabel(r'$\Delta\mathrm{V}$'+" ("+ r"$\mathrm{cm}^3$" +" "+ r"$\mathrm{mol}^{-1})$",fontsize=14)
ax[0].set_xlabel("P (GPa)",fontsize=14)
ax[1].set_xlabel("P (GPa)",fontsize=14)
ax[0].set_ylim([0,12])
ax[1].set_ylim([0,12])
ax[0].set_xlim([0,140])
ax[1].set_xlim([0,140])
ax[0].text(5,10.6,'a',fontsize=16, fontweight='bold')
ax[1].text(5,10.6,'b',fontsize=16, fontweight='bold')
ax[1].tick_params(labelleft=False)

plot_tool.set_major_axis_font(ax[0],14)
plot_tool.set_major_axis_font(ax[1],14)
plot_tool.show_minor_ticks(ax)


if False:
    fig.savefig('out/fig1.pdf',bbox_inches='tight')
