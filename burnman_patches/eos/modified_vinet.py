# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.
"""
added by Jie Deng
"""

import scipy.optimize as opt
import scipy.integrate as integ
from . import equation_of_state as eos
import warnings
from math import exp
import numpy as np

def bulk_modulus(volume, params):
    """
    compute the bulk modulus as per the third order
    Vinet equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in :math:`[Pa]`.
    """

    x = volume / params['V_0']
    eta = (3. / 2.) * (params['Kprime_0'] - 1.)

    K = (params['K_0'] * pow(x, -2. / 3.)) * \
        (1 + ((eta * pow(x, 1. / 3.) + 1.) * (1. - pow(x, 1. / 3.)))) * \
        np.exp(eta * (1. - pow(x, 1. / 3.)))
    return K


def vinet_static(x, params):
    """
    equation for the third order Vinet equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    v  -> p
    x is volume/volume0
    """
    eta         = (3. / 2.) * (params['Kprime_0'] - 1.)
    p_static    = 3. * params['K_0'] * (pow(x, -2. / 3.)) * (1. - (pow(x, 1. / 3.))) \
        * np.exp(eta * (1. - pow(x, 1. / 3.)))
    return p_static
    

def volume(pressure, t, params):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    p, t - > v
    """
    ### This corner case process method is borrowed from pytheos ###
    if pressure <= 1.e-5:
        return params['V_0']
    func = lambda x: vinet_static(x / params['V_0'], params) - pressure
    vp = opt.brentq(func, 0.1 * params['V_0'], 1.5 * params['V_0'])
    ############# Below is different from vinet eos#############
    #############       ref. Komabayashi 2014      #############
    eta_k = vp/ params['V_0']
    alpha = params['alpha0'] * np.exp(-params['delta0']/params['kappa'] * (1-eta_k**params['kappa']))
#    print("vp is",vp)
#    print("alpha is",alpha)
    vt    = vp*np.exp(alpha*(t - 298));
    return vt

def volume_vector(pressure, t, params):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    p, t - > v
    """
    volume_v = np.vectorize(volume, excluded=[1, 2])
    return volume_v(pressure, t, params)

def vinet(x, t, params):
    """
    equation for the third order Vinet equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    v, t  -> p
    here x is volume/volume0
    """
    v = x * params['V_0']
    func = lambda pressure: volume(pressure,t, params) - v
    guess = vinet_static(x, params)
    print("guess is",guess)
    min_guess = guess
    if guess < 0:
        min_guess = 1e5
    p = opt.brentq(func, min_guess, 10 * guess)
    return p

def vinet_vector(x, t, params):
    """
    equation for the third order Vinet equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    v, t  -> p
    here x is volume/volume0
    """
    vinet_v = np.vectorize(vinet, excluded=[1, 2])
    p = vinet_v(x,t,params)
    return p

class MVinet(eos.EquationOfState):

    """
    Base class for the isothermal Vinet equation of state.  This is third order in strain, and
    has no temperature dependence.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
#        print("volume")
        return volume(pressure, temperature, params)

    def pressure(self, temperature, volume, params):
        print("here")
        return vinet_vector(volume / params['V_0'],temperature, params)

    def gibbs_free_energy_slow(self, pressure, temperature, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        This is not consisitent with the the interface gibbs_free_energy defined in eos
        which requires four inputs (pressure, temperature, volume, params):
        ----------
        This cannot be invoked by evaluate function!
        The correct format is as follows,
        fe_fcc.method.gibbs_free_energy(1e9, 3000, fe_fcc.params)

        """
        G1bar = params['a'] + params['b']*temperature + params['c']*temperature*np.log(temperature) + \
                params['d']*(temperature**2) + params['e']/temperature + params['f']*(temperature**0.5);
        if pressure < 1e9:
            bins = 100
        else:
            bins = pressure/1e9*1e2
        pin = np.linspace(1e5,pressure,bins)
        vin = volume_vector(pin, temperature, params)
        G   = np.trapz(vin,pin) + G1bar
        return G
    
    def gibbs_free_energy(self, pressure, temperature, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        This is not consisitent with the the interface gibbs_free_energy defined in eos
        which requires four inputs (pressure, temperature, volume, params)
        pressure and temperature are must be scalar
        example
        ----------
        mineral.method.gibbs_free_energy(2e9, 3000, params)

        """
        G1bar = params['a'] + params['b']*temperature + params['c']*temperature*np.log(temperature) + \
                params['d']*(temperature**2) + params['e']/temperature + params['f']*(temperature**0.5);
        dG = integ.quad(volume_vector,1e5,pressure, args=(temperature,params))
        G  = dG + G1bar
        return G[0]
    
    def gibbs_free_energy_vector_slow(self, pressure, temperature, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        This is not consisitent with the the interface gibbs_free_energy defined in eos
        which requires four inputs (pressure, temperature, volume, params):
        ----------
        This cannot be invoked by evaluate function!
        vectorized gibbs_free_energy
        fe_fcc.method.gibbs_free_energy([1e9, 2e9], 3000, fe_fcc.params)
        ### this is way too slow ####
        ### accelerate it using dynamic programming??? ####
        """
        gout = np.zeros(pressure.shape)
        for i in range(len(pressure)):
            gout[i] = self.gibbs_free_energy_slow(pressure[i], temperature, params)
        return gout

    
    def gibbs_free_energy_vector(self, pressure, temperature, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        This is not consisitent with the the interface gibbs_free_energy defined in eos
        which requires four inputs (pressure, temperature, volume, params):
        pressure can be vector while temperature  must be scalar
        ----------
        This cannot be invoked by evaluate function!
        vectorized gibbs_free_energy
        fe_fcc.method.gibbs_free_energy([1e9, 2e9], 3000, fe_fcc.params)

        """
        gibbs_v = np.vectorize(self.gibbs_free_energy, excluded=[1, 2])
        return gibbs_v(pressure, temperature, params)

    
    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus(volume, params)

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus :math:`K_s` of the mineral. :math:`[Pa]`.
        """
        return bulk_modulus(volume, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        Currently not included in the Vinet EOS, so omitted.
        """
        return 0.

    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[1/K]`
        """
        return 0.

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[unitless]`
        """
        return 0.

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # G is not included in the Vinet EOS so we shall set them to NaN's
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_0' not in params:
            params['Gprime_0'] = float('nan')

        # check that all the required keys are in the dictionary
        expected_keys = ['V_0', 'K_0', 'Kprime_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # now check that the values are reasonable.  I mostly just
        # made up these values from experience, and we are only
        # raising a warning.  Better way to do this? [IR]
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-3:
            warnings.warn('Unusual value for V_0', stacklevel=2)
        if params['K_0'] < 1.e9 or params['K_0'] > 1.e13:
            warnings.warn('Unusual value for K_0', stacklevel=2)
        if params['Kprime_0'] < -5. or params['Kprime_0'] > 10.:
            warnings.warn('Unusual value for Kprime_0', stacklevel=2)
