"""
A class LensInfo that holds relevant information about a
gravitational lens systems and calculate some relevant parameters
"""

import numpy as np
from astropy import constants as const
from math import pi

class LensInfo:

   def __init__(self, nlens):
       """
       Just sets up an empty container for holding lens information

       Input:
         nlens - number of lenses for which information will be stored
       """

       lenstype = [('zl',float),('zs',float),('theta_E',float),('dl',float),\
                       ('ds',float),('dls',float),('r_E',float), \
                       ('log_M_E',float)]
       self.data = np.ones((nlens),dtype=lenstype)
       for i in ('zl','zs','theta_E','dl','ds','dls','r_E','log_M_E'):
           self.data[i] *= -1.

   def calc_M_ein(self, h=0.7):
       """
       Calculates the Einstein ring mass based on an Einstein ring angular
       radius and angular diameter distances for a system
       """

       """ Simplify the constants """
       c = const.c.cgs.value
       G = const.G.cgs.value
       Mpc = 1000. * const.kpc.cgs.value
       Msun = const.M_sun.cgs.value

       """ Conversion factors """
       deg2rad = pi / 180.
       arcsec2rad = deg2rad / 3600.

       """ Mask the points that do not have the required information """
       data = self.data
       mask = (data['zl']>0) & (data['zs']>0) & (data['theta_E']>0)
       D = data['dl'] * data['ds'] / data['dls']

       """ Calculate the Einstein ring mass """
       data['log_M_E'][mask] = \
           np.log10((c * data['theta_E'][mask])**2 * arcsec2rad**2 * D[mask] * Mpc \
                        / (4. * G * Msun * h))

   def calc_R_ein(self):
       """
       Uses the angular diameter distance to the lens, D_l (assumed to have units
       of Mpc), and the angular Einstein ring radius (assumed to be in arcsec)
       to calculate the physical Einstein ring radius in kpc
       """

       """ Set the mask for valid data """
       data = self.data
       mask = (data['dl']>0) & (data['theta_E']>0)

       """ 
       Do the conversion 
       The factor of 1000. is to convert the units from Mpc to kpc
       """
       deg2rad = pi / 180.
       data['r_E'][mask] = 1000. * data['dl'][mask] * data['theta_E'][mask] * \
           deg2rad / 3600.
