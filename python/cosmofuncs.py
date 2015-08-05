"""
Functions for calculating cosmological parameters, all contained in
the Cosmodat class.

Very much under construction.
"""

import numpy as n
import os

class Cosmodat:

   def __init__(self, filename='interactive', zcol=0):
      """
      Reads in redshift(s) for which cosmology will be calculated.

      Inputs:
         filename:  Input filename containing a list of redshifts (one per
                    line).  The default, 'interactive' is for interactive
                    input of a single redshift.
         zcol:      Column in the input file (if any) containing the redshifts
                    Default is zcol=0 (first column)
      """

      if filename != 'interactive':
         try:
            self.z = n.loadtxt(filename,unpack=True,usecols=(zcol,))
         except:
            print ''
            print 'ERROR: File %s was not found or data could not be read' % \
                filename
            return
      else:
         z = float(raw_input('Enter redshift: '))
         """ Convert to a numpy array """
         self.z = n.atleast_1d(z)

      self.cosdat = None

   #-----------------------------------------------------------------------

   def get_cosmo(self):
      """
      Interactively loads a cosmology into the class.

      """

      self.omega_m  = omega_m
      self.omega_de = omega_de
      self.w        = w

   #-----------------------------------------------------------------------

   def calc_cosdist(self, zmin=1.0e-30,zmax=100.,new_cosmo=False, verbose=True):
      """
      Calculates cosmological distances, given the input cosmological
      world model.

      *** For now, does this by calling the C program cosmocalc.c ***
      *** This means that for now the input cosmology is ignored since
          these parameters are set interactively by cosmocalc.c
      ***

      Inputs:
         omega_m:   Omega_matter
         omega_de:  Omega_dark-energy
         w:         Dark energy equation of state parameter
      """

      """ Start by writing out a temporary file containing the redshifts """
      self.zgood = self.z[(self.z>=zmin) & (self.z<=zmax)]
      n.savetxt('tmp_z.txt',self.zgood)

      """
      Get cosmological distances by running cosmocalc.c
      """
      print ""
      print "Getting cosmological info by running cosmocalc.c"
      print "------------------------------------------------"
      os.system('cosmocalc -b tmp_z.txt tmp_cosmo.out')
      cosdat = n.atleast_2d(n.loadtxt('tmp_cosmo.out'))
      self.omega_m  = cosdat[:,1]
      self.omega_de = cosdat[:,2]
      self.w        = cosdat[:,3]
      self.hz       = cosdat[:,4]
      self.tlook    = cosdat[:,5]
      self.DM       = cosdat[:,6]
      self.D_A      = cosdat[:,7]
      self.D_L      = cosdat[:,8]
      self.D_L_area = cosdat[:,9]
      self.asec     = cosdat[:,10]
      self.Mpc      = cosdat[:,11]
      self.Mpc_com  = cosdat[:,12]

      if verbose:
         print "  z   Om_m Om_L  w   H(z)    t_l     DM   D_A  D_L  4*PI*DL^2 1asec  1Mpc  1Mpc"
         print "----- ---- ---- ---- ---- --------- ----- ---- ---- --------- ----- ----- -----"
         for i in range(self.zgood.size):
            print "%5.3f %4.2f %4.2f %4.1f %4.0f %9.3e %5.2f %4.0f %4.0f %9.3e %5.2f %5.2f %5.2f" \
                %(self.zgood[i],self.omega_m[i],self.omega_de[i],self.w[i],
                  self.hz[i],self.tlook[i],self.DM[i],self.D_A[i],self.D_L[i],
                  self.D_L_area[i],self.asec[i],self.Mpc[i],self.Mpc_com[i])


