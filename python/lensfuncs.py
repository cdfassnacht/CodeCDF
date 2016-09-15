"""
A class LensInfo that holds relevant information about a
gravitational lens systems and calculate some relevant parameters
"""

import numpy as np
from astropy import constants as const
from math import pi
from matplotlib import pyplot as plt

class LensInfo:

   def __init__(self, intable, h=0.7):
       """
       Sets up the LensInfo class and populates it with the data in the
        input table.
       Note: the input table is expected to be an astropy.table.Table class,
        that may very well have been created by a call to astropy.ascii.read

       Input:
         intable - input table, 
       """

       self.data = intable

       """ Create missing columns if needed """
       for i in ('zl','zs','thetaE','D_l','D_s','D_ls','R_E','logM_E','m_814',
                 'logM_proj'):
          try:
             tmp = self.data[i]
          except:
             self.data[i] = np.ones(len(self.data)) * -1.

       """ 
       Set up masks.  The logM_E and R_E masks are set to None here, but
        can be set later either by the relevant calc_* function, or by
        a call outside the class
       """
       self.zlmask = self.data['zl']>0
       self.zsmask = self.data['zs']>0
       self.remask = None
       self.memask = None

       """ Set plotting parameters to default values """
       self.label = None
       self.color = 'b'
       self.alpha = 1.0
       self.symb  = 'bo'
       self.ls    = 'solid'
       self.opensymb = False

   #---------------------------------------------------------------------------

   def plot_hist(self, param, mask=None, bins=None, xlab=None, pdf=False):
      """

      Generic routine for plotting a histogram that replaces all of the 
       parameter-specific routines that were used before.

      NOTE: This routine uses the color, alpha, opensymb, and ls parameters
       that are stored within the class rather than setting them on the command
       line.  The default values for these parameters are set in the __init__
       method.  If other values are desired, then they need to be changed
       explicitly BEFORE calling plot_hist.
       For example, to create an instance of the class and then reset, e.g.,
        the color from its default value before calling plot_hist,
        you can do something like the following:

           mysamp = LensInfo(intable)
           mysamp.color = 'r'
           mysamp.plot_hist('zl')

      """


      """ Set the data """
      if mask is None:
         data = self.data[param]
      else:
         data = self.data[param][mask]

      """ Set the histogram type """
      if self.opensymb:
         histtype = 'step'
         lw = 3
      else:
         histtype = 'bar'
         lw = None

      """ Plot the histogram """
      if bins is None:
         tmp = plt.hist(data,color=self.color,alpha=self.alpha,histtype=histtype,
                        lw=lw,ls=self.ls,label=self.label)
         self.bins = tmp[1]
      else:
         self.bins = bins
         tmp = plt.hist(data,bins=bins,color=self.color,alpha=self.alpha,
                        histtype=histtype,lw=lw,ls=self.ls,label=self.label)
      self.hist = tmp[0]
      if xlab:
         plt.xlabel(xlab)
      plt.ylabel('Relative numbers')

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
      plt.plot(data['zl'][mask],data['logM_E'][mask],symb,ms=size,
               fillstyle=fs,alpha=alpha)
      plt.xlabel('Lens redshift')
      plt.ylabel(r'Projected mass inside $R_{\rm Ein}$ (h=%3.1f)' % h)

   #---------------------------------------------------------------------------

   def plot_M_proj(self, Rproj, h=0.7, symb='bo', size=5, alpha=1.,
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
      plt.plot(self.data['zl'][self.memask],self.logM_proj,symb,ms=size,
               fillstyle=fs,alpha=alpha)
      plt.xlabel('Lens redshift')
      plt.ylabel('Projected mass inside R = %4.1f kpc (h=%3.1f)' % (Rproj,h))

   #---------------------------------------------------------------------------

   def hist_M_ein(self, h=0.7, color='b', alpha=1., opensymb=False,
                  label=None, ls='solid'):
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
      plt.hist(data['logM_E'][mask],color=color,label=label)

   #---------------------------------------------------------------------------

   def hist_zlens(self, color='b', alpha=1., opensymb=False,
                  label=None, ls='solid', bins=None):

      """ Now make the plot """
      if opensymb:
         histtype = 'step'
         lw = 3
      else:
         histtype = 'bar'
         lw = None
      #print histtype
      if bins is None:
         tmp = plt.hist(self.data['zl'][self.zlmask],color=color,alpha=alpha,
                        histtype=histtype,lw=lw,ls=ls,label=label)
         self.bins = tmp[1]
      else:
         self.bins = bins
         plt.hist(self.data['zl'][self.zlmask],bins=bins,color=color,alpha=alpha,
                  histtype=histtype,lw=lw,ls=ls,label=label)

   #---------------------------------------------------------------------------

   def hist_R_E(self, color='b', alpha=1., opensymb=False,
                label=None, ls='solid', bins=None):

      """ Now make the plot """
      if opensymb:
         histtype = 'step'
         lw = 3
      else:
         histtype = 'bar'
         lw = None
      #print histtype
      if bins is None:
         tmp = plt.hist(self.data['R_E'][self.remask],color=color,alpha=alpha,
                        histtype=histtype,lw=lw,ls=ls,label=label)
         self.bins = tmp[1]
      else:
         self.bins = bins
         plt.hist(self.data['M_E'][self.remask],bins=bins,color=color,
                  alpha=alpha,histtype=histtype,lw=lw,ls=ls,label=label)

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
      plt.hist(self.logM_proj,color=color,alpha=alpha,histtype=histtype,
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

   def calc_M_proj(self, Rproj, h=0.7):
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
      logMproj = log_scale + self.data['logM_E'][self.memask]
      self.data['logM_proj'][self.memask] = logMproj
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
       self.memask = (data['zl']>0) & (data['zs']>0) & (data['thetaE']>0)
       D = data['D_l'] * data['D_s'] / data['D_ls']

       """ Combine the constants into one factor """
       cfac = (c * arcsec2rad)**2 * Mpc / (4. * G * Msun * h)
       """ Calculate the Einstein ring mass """
       data['logM_E'][self.memask] = \
           np.log10(cfac * (data['thetaE'][self.memask])**2 * D[self.memask])
                        

   #---------------------------------------------------------------------------

   def calc_R_ein(self, h=0.7):
       """
       Uses the angular diameter distance to the lens, D_l (assumed to have units
       of Mpc), and the angular Einstein ring radius (assumed to be in arcsec)
       to calculate the physical Einstein ring radius in kpc
       """

       """ Set the mask for valid data """
       data = self.data
       mask = (data['D_l']>0) & (data['thetaE']>0)

       """ 
       Do the conversion 
       The factor of 1000. is to convert the units from Mpc to kpc
       """
       deg2rad = pi / 180.
       data['R_E'][mask] = 1000. * data['D_l'][mask] * data['thetaE'][mask] * \
           deg2rad / (3600. * h)
       self.remask = mask
