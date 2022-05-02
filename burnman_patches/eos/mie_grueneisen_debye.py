# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
import scipy.optimize as opt
import warnings

from . import equation_of_state as eos
from . import birch_murnaghan as bm
from . import debye
from .. import constants
from ..tools import bracket


class MGDBase(eos.EquationOfState):

    """
    Base class for a generic finite-strain Mie-Grueneisen-Debye
    equation of state.  References for this can be found in many
    places, such as Shim, Duffy and Kenichi (2002) and Jackson and Rigden
    (1996).  Here we mostly follow the appendices of Matas et al (2007).
    Of particular note is the thermal correction to the shear modulus, which
    was developed by Hama and Suito (1998).
    """

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume (EQ B6)
        """
        return self._grueneisen_parameter(params['V_0'] / volume, params)

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        EQ B7
        """
        T_0 = params['T_0']
        func = lambda x: bm.birch_murnaghan(params['V_0'] / x, params) + \
            self._thermal_pressure(temperature, x, params) - \
            self._thermal_pressure(T_0, x, params) - pressure
        try:
            sol = bracket(func, params['V_0'], 1.e-2 * params['V_0'])
        except:
            raise ValueError(
                'Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?')
        return opt.brentq(func, sol[0], sol[1])

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ B8
        """
        T_0 = params['T_0']
        K_T = bm.bulk_modulus(volume, params) + \
            self._thermal_bulk_modulus(temperature, volume, params) - \
            self._thermal_bulk_modulus(T_0, volume, params)  # EQB13
        return K_T

    # calculate the mgd shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ B11
        """
        T_0 = params['T_0']
        if self.order == 2:
            return bm.shear_modulus_second_order(volume, params) + \
                self._thermal_shear_modulus(temperature, volume, params) - \
                self._thermal_shear_modulus(T_0, volume, params)  # EQ B11
        elif self.order == 3:
            return bm.shear_modulus_third_order(volume, params) + \
                self._thermal_shear_modulus(temperature, volume, params) - \
                self._thermal_shear_modulus(T_0, volume, params)  # EQ B11
        else:
            raise NotImplementedError("")

    # heat capacity at constant volume
    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature, and volume [J/K/mol]
        """
        Debye_T = self._debye_temperature(params['V_0'] / volume, params)
        C_v = debye.heat_capacity_v(temperature, Debye_T, params['n'])
        return C_v

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        """
        C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        gr = self._grueneisen_parameter(params['V_0'] / volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        alpha = gr * C_v / K / volume
        return alpha

    # heat capacity at constant pressure
    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self._grueneisen_parameter(params['V_0'] / volume, params)
        C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v * (1. + gr * alpha * temperature)
        return C_p

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ D6
        """
        K_T = self.isothermal_bulk_modulus(
            pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self._grueneisen_parameter(params['V_0'] / volume, params)
        K_S = K_T * (1. + gr * alpha * temperature)
        return K_S

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        EQ B7
        """
        T_0 = params['T_0']
        return bm.birch_murnaghan(params['V_0'] / volume, params) + \
            self._thermal_pressure(temperature, volume, params) - \
            self._thermal_pressure(T_0, volume, params)

    # calculate the thermal correction to the shear modulus as a function of
    # V, T
    def _thermal_shear_modulus(self, T, V, params):
        if T > 1.e-10:
            gr = self._grueneisen_parameter(params['V_0'] / V, params)
            Debye_T = self._debye_temperature(params['V_0'] / V, params)
            G_th = 3. / 5. * (self._thermal_bulk_modulus(T, V, params) -
                              6 * constants.gas_constant * T * params['n'] / V * gr * debye.debye_fn(Debye_T / T))  # EQ B10
            return G_th
        else:
            return 0.

    # compute the Debye temperature in K.  Takes the
    # parameter x, which is V_0/V (molar volumes).
    # Depends on the reference grueneisen parameter,
    # the reference Debye temperature, and the factor
    # q_0, see Matas eq B6
    def _debye_temperature(self, x, params):
        return params['Debye_0'] * np.exp((params['grueneisen_0'] -
                                           self._grueneisen_parameter(x, params)) / params['q_0'])

    # compute the grueneisen parameter with depth, according
    # to q_0.  Takes x=V_0/V. See Matas eq B6
    def _grueneisen_parameter(self, x, params):
#        print("_grueneisen_parameter is", params['grueneisen_0'] * pow(1. / x, params['q_0']))
        return params['grueneisen_0'] * pow(1. / x, params['q_0'])

    # calculate isotropic thermal pressure, see
    # Matas et. al. (2007) eq B4
    def _thermal_pressure(self, T, V, params):
        Debye_T = self._debye_temperature(params['V_0'] / V, params)
        gr = self._grueneisen_parameter(params['V_0'] / V, params)
        P_th = gr * debye.thermal_energy(T, Debye_T, params['n']) / V
        return P_th

    # calculate the thermal correction for the mgd
    # bulk modulus (see matas et al, 2007)
    def _thermal_bulk_modulus(self, T, V, params):
        if T > 1.e-10:
            gr = self._grueneisen_parameter(params['V_0'] / V, params)
            Debye_T = self._debye_temperature(params['V_0'] / V, params)
            K_th = 3. * params['n'] * constants.gas_constant * T / V * gr * \
                ((1. - params['q_0'] - 3. * gr) * debye.debye_fn(
                 Debye_T / T) + 3. * gr * (Debye_T / T) / (np.exp(Debye_T / T) - 1.))  # EQ B5
            return K_th
        else:
            return 0.
        
    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if 'T_0' not in params:
            params['T_0'] = 300.
        
        # First, let's check the EoS parameters for Tref
        bm.BirchMurnaghanBase.validate_parameters(
            bm.BirchMurnaghanBase(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ['molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Finally, check that the values are reasonable.
        if params['T_0'] < 0.:
            warnings.warn('Unusual value for T_0', stacklevel=2)
        if params['molar_mass'] < 0.001 or params['molar_mass'] > 1.:
            warnings.warn('Unusual value for molar_mass', stacklevel=2)
        if params['n'] < 1. or params['n'] > 1000.:
            warnings.warn('Unusual value for n', stacklevel=2)
        if params['Debye_0'] < 1. or params['Debye_0'] > 10000.:
            warnings.warn('Unusual value for Debye_0', stacklevel=2)
        if params['grueneisen_0'] < 0. or params['grueneisen_0'] > 10.:
            warnings.warn('Unusual value for grueneisen_0', stacklevel=2)
        if params['q_0'] < -10. or params['q_0'] > 10.:
            warnings.warn('Unusual value for q_0', stacklevel=2)

        ########### add by Jie to complete MGD eos #########
    def _helmholtz_free_energy_thermal(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        ref. Liebske & Frost 2012 EPSL
        """
        if not 'F_0'  in params:
            Fth = float('nan')
            
        if 'F_0' in params:
        ## solid   eqS7
            x = params['V_0'] / volume
            debyeT = params['Debye_0']*np.exp((params['grueneisen_0'] - self._grueneisen_parameter(x, params))/params['q_0'])
#            print("debyeT is",debyeT)
            Fth = debye.helmholtz_free_energy(temperature, debyeT, params['n'])
        return Fth    
        
    def _helmholtz_free_energy_cold(self, pressure, temperature, volume, params):
        x = params['V_0'] / volume
        a = 3/2*(params['Kprime_0'] -4)
        f = 1. / 2. * (pow(x, 2. / 3.) - 1.)
        Fcold = 9*params['K_0']*params['V_0']*(1/2*(f**2)+1/3*a*(f**3))
#        print("F_cold is", Fcold)
        return Fcold

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        if not 'F_0' in params:
            raise NotImplementedError("") 
        Fcold = self._helmholtz_free_energy_cold(pressure, temperature, volume, params)
        Fth_T = self._helmholtz_free_energy_thermal(pressure, temperature, volume, params)
        Fth_T0 = self._helmholtz_free_energy_thermal(pressure, params['T_0'], volume, params)
        return params['F_0'] + Fcold + Fth_T - Fth_T0
                
    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        """
        G = self.helmholtz_free_energy(
            pressure, temperature, volume, params) + pressure * volume
        return G

    def internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy at the pressure and temperature of the mineral [J/mol]
        """
        return self.helmholtz_free_energy( pressure, temperature, volume, params) + \
            temperature * \
            self.entropy(pressure, temperature, volume, params)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """
        x = params['V_0'] / volume
        f = 1. / 2. * (pow(x, 2. / 3.) - 1.)
        Debye_T = self._debye_temperature(params['V_0'] / volume, params)
        S = debye.entropy(temperature, Debye_T, params['n'])
        return S

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature of the mineral [J/mol]
        """

        return self.helmholtz_free_energy( pressure, temperature, volume, params) + \
            temperature * self.entropy( pressure, temperature, volume, params) + \
            pressure * volume
        ########### add by Jie to complete MGD eos #########


class MGD3(MGDBase):

    """
    MGD equation of state with third order finite strain expansion for the
    shear modulus (this should be preferred, as it is more thermodynamically
    consistent.
    """

    def __init__(self):
        self.order = 3


class MGD2(MGDBase):

    """
    MGD equation of state with second order finite strain expansion for the
    shear modulus.  In general, this should not be used, but sometimes
    shear modulus data is fit to a second order equation of state.  In that
    case, you should use this.  The moral is, be careful!
    """

    def __init__(self):
        self.order = 2

########## liquid eos based on Liebske and Frost, 2012 EPSL ##########
##########                      add by Jie                  ##########
class LF12L(MGDBase):
    """
    works only for liquid
    """
    def __init__(self):
        self.order = 3
        
    def _grueneisen_parameter(self, x, params):
        gru = params['grueneisen_0'] + params['grueneisen_prime'] * (1. / x -1)
#        print("grueneisen_parameter is", gru)
        return gru
    
    def _thermal_pressure(self, temperature, volume, params):
        # liquid eq S5b
        P_th = self._grueneisen_parameter(params['V_0'] / volume, params)/volume * \
                params['Cv']*(temperature - params['T_0'])
        return P_th

    def _helmholtz_free_energy_thermal(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        ref. Liebske & Frost 2012 EPSL
        """
        if not 'S_0'  in params:
            Fth = float('nan')
            
        if 'S_0' in params:
        ## liquid   eqS7
            dT = temperature - params['T_0']
            x  = params['V_0'] / volume
            dg = self._grueneisen_parameter(x,params) - params['grueneisen_0']
            Fth = -params['S_0'] * dT - params['Cv'] * (temperature*np.log(temperature/params['T_0']) - dT) - \
                  params['Cv'] * dT * ((params['grueneisen_0'] -  params['grueneisen_prime'])*np.log(volume/params['V_0']) + dg)
#        print("Fth is", Fth)
        return Fth

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        return params['F_0'] + self._helmholtz_free_energy_cold(pressure, temperature, volume, params) + \
                self._helmholtz_free_energy_thermal(pressure, temperature, volume, params)
    
    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if 'T_0' not in params:
            params['T_0'] = 300.
        
        # First, let's check the EoS parameters for Tref
        bm.BirchMurnaghanBase.validate_parameters(
            bm.BirchMurnaghanBase(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
#        expected_keys = ['molar_mass', 'n', 'Debye_0', 'grueneisen_0']# q_0 not expected for liquid added by Jie
#        for k in expected_keys:
#            if k not in params:
#                raise KeyError('params object missing parameter : ' + k)

        # Finally, check that the values are reasonable.
#        if params['T_0'] < 0.:
#            warnings.warn('Unusual value for T_0', stacklevel=2)
#        if params['molar_mass'] < 0.001 or params['molar_mass'] > 1.:
#            warnings.warn('Unusual value for molar_mass', stacklevel=2)
#        if params['n'] < 1. or params['n'] > 1000.:
#            warnings.warn('Unusual value for n', stacklevel=2)
#        if params['Debye_0'] < 1. or params['Debye_0'] > 10000.:
#            warnings.warn('Unusual value for Debye_0', stacklevel=2)
#        if params['grueneisen_0'] < 0. or params['grueneisen_0'] > 10.:
#            warnings.warn('Unusual value for grueneisen_0', stacklevel=2)

        