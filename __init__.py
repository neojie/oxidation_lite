"""
oxidation
======
source code to generate figures of 
Deng, J., Du, Z., Ghosh, D, Karki, B. & Lee, K.K.M., 2020. \
A magma ocean origin of the divergent redox evolutions of \
rocky planetary bodies and early atmospheres. Nature Communications, 11(1), 2007 

"""
"""
Pre-requisite
======
pandas
lmfit
numpy
scipy : interp1d, opt
burnman : modified `burnman` with `LF_2012` 
uncertainties
matplotlib

"""
"""
Code List
======
/src
fig1   : dV, eos
         tool => cal_dV_this_study
         oxidation_5.xlsx
         
fig2   : redox vs. depth of MO
         tool => cal_dV_this_study
         tool2 => X_cal_vector, fo2_cal_vector

fig3   : redox vs. depth of MO base
         tool => cal_dV_this_study
         tool2 => X_cal_vector, fo2_cal_vector

## codes to generate the supplementary figures are available upon request

---------
A recent paper by borisov et al., 2018 presented some more data at 1 bar, which 
we may use to construct better model for activity ratio. Also, we can gest our 
best-fitting model against these new data ( 195 new data).

ref: Ferric/ferrous ratio in silicate melts: a new model for 1 atm data

"""

#### import src folder 
from . import src
from .src import *


