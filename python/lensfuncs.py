"""
A class LensInfo that holds relevant information about a
gravitational lens systems and calculate some relevant parameters
"""

import numpy as np
from astropy import constants as const
from math import pi
from matplotlib import pyplot as plt

class LensInfo:

   def __init__(self, nlens):
       """
       Just sets up an empty container for holding lens information on
        a set of nlens lens systems

       Input:
         nlens - number of lenses for which information will be stored
       """

       lenstype = [('name','S20'),('zl',float),('zs',float),('theta_E',float), \
                       ('dl',float),('ds',float),('dls',float), \
                       ('R_E',float),('log_M_E',float)]
       self.data = np.ones((nlens),dtype=lenstype)
       for i in ('zl','zs','theta_E','dl','ds','dls','R_E','log_M_E'):
           self.data[i] *= -1.

       self.memask = None
       self.remask = None

   #---------------------------------------------------------------------------

   def plot_M_ein(self, h=0.7, symb='bo',size=5,alpha=1.,opensymb=False):
      """
      Makes a plot of Einstein ring radius vs. lens redshift for systems
      that have the appropriate information.
      """

      """ First calculate the Einstein ring masses if needed """
      if self.memask is None:
         self.calc_M_ein(h)

      """ Now make the plot """
      data = self.data
      mask = self.memask
      if opensymb:
         fs = 'none'
      else:
         fs = 'full'
      plt.plot(data['zl'][mask],data['log_M_E'][mask],symb,ms=size,
               fillstyle=fs,alpha=alpha)
      plt.xlabel('Lens redshift')
      plt.ylabel(r'Projected mass inside $R_{\rm Ein}$ (h=%3.1f)' % h)

   #---------------------------------------------------------------------------

   def plot_M_proj(self, Rproj, h=0.7, symb='o', color='b', size=5, alpha=1.,
                   opensymb=False, label=None):
      """
      Makes a plot of projeced mass inside a ring of a given radius vs. lens 
       redshift for systems that have the appropriate information.
      The masses are calculated by first calculating the Einstein ring projected
       masses and also the Einstein radii in h^-1 kpc.  Then assume an isothermal
       profile, which means that projected mass within a circle centered on
       the lens scales as Rproj, to extrapolate to the mass within the given
       value of Rproj
      """

      """ Calculate the projected masses within the requested radius """
      self.calc_R_proj(Rproj,h)

      """ Now make the plot """
      if opensymb:
         fs = 'none'
      else:
         fs = 'full'
      plt.plot(self.data['zl'][self.memask],self.log_M_proj,symb,color=color,
               ms=size,fillstyle=fs,alpha=alpha)
      plt.xlabel('Lens redshift')
      plt.ylabel('Projected mass inside R = %4.1f kpc (h=%3.1f)' % (Rproj,h))

   #---------------------------------------------------------------------------

   def hist_M_ein(self, h=0.7, color='b', alpha=1., opensymb=False,
                  label=None):
      """ First calculate the Einstein ring masses if needed """
      if self.memask is None:
         self.calc_M_ein(h)

      """ Now make the plot """
      data = self.data
      mask = self.memask
      if opensymb:
         fs = 'none'
      else:
         fs = 'full'
      plt.hist(data['log_M_E'][mask],color=color,label=label)

   #---------------------------------------------------------------------------

   def hist_M_proj(self, Rproj, h=0.7, color='b', alpha=1., opensymb=False,
                   label=None, ls='solid'):

      """ Calculate the projected masses within the requested radius """
      self.calc_R_proj(Rproj,h)

      """ Now make the plot """
      if opensymb:
         histtype = 'step'
         lw = 3
      else:
         histtype = 'bar'
         lw = None
      #print histtype
      plt.hist(self.log_M_proj,color=color,alpha=alpha,histtype=histtype,
               lw=lw,ls=ls,label=label)

   #---------------------------------------------------------------------------

   def plot_R_ein(self, h=0.7, symb='bo',size=5,alpha=1.,opensymb=False):
      """
      Makes a plot of Einstein ring radius vs. lens redshift for systems
      that have the appropriate information.
      """

      """ First calculate the Einstein ring radii if needed """
      if self.remask is None:
         self.calc_R_ein(h)

      """ Now make the plot """
      data = self.data
      mask = self.remask
      if opensymb:
         fs = 'none'
      else:
         fs = 'full'
      plt.plot(data['zl'][mask],data['R_E'][mask],symb,ms=size,
               fillstyle=fs,alpha=alpha)
      plt.xlabel('Lens redshift')
      plt.ylabel('Einstein ring radius (kpc, h=%3.1f)' % h)

   #---------------------------------------------------------------------------

   def calc_R_proj(self, Rproj, h=0.7):
      """
      Makes a plot of projeced mass inside a ring of a given radius vs. lens 
       redshift for systems that have the appropriate information.
      The masses are calculated by first calculating the Einstein ring projected
       masses and also the Einstein radii in h^-1 kpc.  Then assume an isothermal
       profile, which means that projected mass within a circle centered on
       the lens scales as Rproj, to extrapolate to the mass within the given
       value of Rproj
      """

      """ First calculate the Einstein ring radii and masses if needed """
      if self.remask is None:
         self.calc_R_ein(h)
      if self.memask is None:
         self.calc_M_ein(h)

      """ 
      Now scale each of the masses to extrapolate to the given value of Rproj.
      Since the mass is stored as a log10, and the scaling is just a ratio
       that will be used to multiply by the mass, we can just add the log
       of the ratio
      """
      log_scale =  np.log10(Rproj / self.data['R_E'][self.memask])
      log_Mproj = log_scale + self.data['log_M_E'][self.memask]
      self.log_M_proj = log_Mproj
      self.R_proj = Rproj

   #---------------------------------------------------------------------------

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
       self.memask = (data['zl']>0) & (data['zs']>0) & (data['theta_E']>0)
       D = data['dl'] * data['ds'] / data['dls']

       """ Combine the constants into one factor """
       cfac = (c * arcsec2rad)**2 * Mpc / (4. * G * Msun * h)
       """ Calculate the Einstein ring mass """
       data['log_M_E'][self.memask] = \
           np.log10(cfac * (data['theta_E'][self.memask])**2 * D[self.memask])
                        

   #---------------------------------------------------------------------------

   def calc_R_ein(self, h=0.7):
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
       data['R_E'][mask] = 1000. * data['dl'][mask] * data['theta_E'][mask] * \
           deg2rad / (3600. * h)
       self.remask = mask
