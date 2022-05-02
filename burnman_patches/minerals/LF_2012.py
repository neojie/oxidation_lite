# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
added by Jie Deng 
Mar 20 2019
^^^^^^^^^^^^^^^^^^

Minerals from Liebske and Frost (2002)


"""
from __future__ import absolute_import
#from .. import mineral_helpers as helpers
#from ..mineral import Mineral
from burnman import Mineral

        
class mgsio3(Mineral):

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


class mg_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd3',
            'T_0': 300,
            'V_0': 24.25e-6,  
            'K_0': 251e9,
            'Kprime_0': 4.1,
            'molar_mass': .1,
            'n': 5, 
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'Debye_0': 905,
            'F_0': -1368000}

        Mineral.__init__(self)

mg_bridgmanite = mg_perovskite
