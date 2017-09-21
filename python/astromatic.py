"""
astromatic.py

A library containing functions that run files to prepare for or to invoke
the cataloging and astrometry programs provided by astromatic.net (formerly
terapix).  

Instrument-specific functions
-----------------------------
 make_cat_irac
 make_cat_wfc3
 make_cat_fors2
 make_cat_isaac
 make_cat_niri
 make_cat_suprimecam
 make_cat_hawki
 make_cat_vircam (use for VISTA VHS)

Generic functions
------------------
 make_fits_cat   - the workhorse catalog-generating function (uses SExtractor)
 make_reg_file   - converts a catalog into a ds9 region file
 run_scamp       - runs the Astromatic scamp code to solve for the WCS
 run_swarp       - runs the Astromatic swarp code to resample onto a new WCS
 do_photom       - NOT YET IMPLEMENTED.

"""

import os
import numpy as n
from astropy.io import fits as pf
import catfuncs as cf
import imfuncs as imf

# =======================================================================

class fits2cat(imf.Image):
    """
    """

    def __init__(self, fitsfile, gain='gain', texp='exptime', hext=0, 
                 verbose=True):
        """
        """

        """ Link to the superclass """
        self.fitsfile = fitsfile
        imf.Image.__init__(self, fitsfile, verbose=verbose)

        """ Set up default values """
        self.outcat = fitsfile.replace('.fits', '.cat')
        self.regfile = fitsfile.replace('.fits', '.reg')
        self.whtfile = fitsfile.replace('.fits', '_wht.fits')
        self.config = 'sext_astfile.config'
        self.gain = None
        self.texp = None
        self.namecol = 0
        self.racol = 1
        self.deccol = 2

        """ Get the gain and exposure time """
        self.set_gain(gain=gain, hext=hext, verbose=verbose)
        self.set_texp(texp=texp, hext=hext, verbose=verbose)

    # -----------------------------------------------------------------------

    def set_gain(self, gain='gain', hext=0, verbose=True):
        """
        Sets the gain.  This can be done in one of two ways:
         1. If the gain parameter is is a string, then the gain will be read
            from the fits header card identified by that string.
         2. If the gain parameter is a number, then the gain will be set to
            that value
        If option 1 is chosen but the request information is not in the
         fits header, then the gain is set to 1.0
        """

        if type(gain) == float:
            self.gain = gain
        elif type(gain) == str:
            hdr = self.hdu[hext].header
            if gain in hdr.keys():
                self.gain = hdr[gain]

        """ If self.gain is still None, then it must not have been found """
        if self.gain is None:
            if verbose:
                print('Warning: gain not found.  Setting gain to 1.0')
            self.gain = 1.0

    # -----------------------------------------------------------------------

    def set_texp(self, texp='exptime', hext=0, verbose=True):
        """
        Sets the exposure time.  This can be done in one of two ways:
         1. If the texp parameter is a string, then the exposure time will be
            read from the fits header card identified by that string.
         2. If the texp paramter is a number, then the gain will be set to that
            value
        If option 1 is chosen but the request information is not in the
         fits header, then the exposure time is set to 1.0
        """

        if type(texp) == float:
            self.texp = texp
        elif type(texp) == str:
            hdr = self.hdu[hext].header
            if texp in hdr.keys():
                self.texp = hdr[texp]

        """ If self.texp is still None, then it must not have been found """
        if self.texp is None:
            if verbose:
                print('Warning: texp not found.  Setting texp to 1.0')
            self.texp = 1.0

    # -----------------------------------------------------------------------

    def make_secat(self, outcat='default', configfile='default',
                   ncoadd=1, satur=64000., zeropt=None,
                   catformat='ldac', det_area=-1, det_thresh=-1., seeing=0.0,
                   whtfile=None, weight_type='MAP_WEIGHT', weight_thresh=None,
                   regfile='default', flag_file=None, logfile=None,
                   verbose=True, racol=None, deccol=None, fluxcol=None,
                   fluxerrcol='fluxerr_auto'):

        """ Set up the configuration file """
        if configfile == 'default':
            configfile = self.config

        """
        Set up the gain and format of the output catalog (default is LDAC)
        """
        if catformat.lower() == 'ascii':
            cattype = 'ASCII_HEAD'
        else:
            cattype = 'FITS_LDAC'
        gain_eff = self.gain * self.texp * ncoadd

        """ Prepare other parameters """
        satur_eff = satur * ncoadd

        """ Prepare the optional portions of the call to SExtractor """
        if outcat == 'default':
            outcat = self.outcat
        sopts = '-GAIN %f -CATALOG_NAME %s -CATALOG_TYPE %s ' \
            % (gain_eff, outcat, cattype)
        sopts += '-SATUR_LEVEL %f ' % satur_eff
        if det_area>0:
            sopts += '-DETECT_MINAREA %d ' % det_area
        if det_thresh>0.:
            sopts += '-DETECT_THRESH %5.2f -ANALYSIS_THRESH %5.2f ' \
                % (det_thresh, det_thresh)
        if zeropt is not None:
            sopts += '-MAG_ZEROPOINT %5.2f ' % zeropt
        if seeing>0.0:
            sopts += '-SEEING_FWHM %6.3f ' % seeing
        if whtfile is None:
            sopts += '-WEIGHT_TYPE NONE '
        else:
            if whtfile == 'default':
                whtfile = self.whtfile
            sopts += '-WEIGHT_TYPE %s -WEIGHT_IMAGE %s ' % \
                (weight_type, whtfile)
        if weight_thresh is not None:
            sopts += '-WEIGHT_THRESH %9.2f ' % weight_thresh
        if flag_file is not None:
            sopts += '-FLAG_IMAGE %s ' % flag_file
        if verbose is False:
            sopts += '-VERBOSE_TYPE QUIET '

        """ Run SExtractor """
        if verbose:
            print('')
            print('Running SExtractor on %s' % self.fitsfile)
            print('Configuration file: %s' % configfile)
            print('Override variables:')
            print sopts
            print('')
        if logfile is None:
            try:
                os.system('sex -c %s %s %s' % (configfile, self.fitsfile, sopts))
            except:
                print ""
                print "ERROR.  Could not run SExtractor on %s" % self.fitsfile
                print ""
                return
        else:
            try:
                os.system('sex -c %s %s %s > %s' % \
                              (configfile, self.fitsfile, sopts, logfile))
            except:
                print ""
                print "ERROR.  Could not run SExtractor on %s" % self.fitsfile
                print ""
                return
        if verbose:
            print ""
            print "Ran SExtractor on %s to produce output catalog %s" % \
                (self.fitsfile, outcat)

        if regfile is not None:
            if catformat=='ascii':
                if racol is not None:
                    self.racol = racol
                if deccol is None:
                    self.deccol = deccol
                tmpcat = cf.Secat(outcat, catformat, racol=self.racol, 
                                  deccol=self.deccol, namecol=self.namecol)
            else:
                tmpcat = cf.Secat(outcat, catformat)

            if regfile == 'default':
                regfile = self.regfile
            tmpcat.make_reg_file(regfile, 1.5, fluxcol=fluxcol,
                                 fluxerrcol=fluxerrcol)
        if verbose:
            print('')

# -----------------------------------------------------------------------

def make_cat_generic(fitsfile, outcat='default', regfile='default',
                     configfile='sext_astfile.config', 
                     whtfile=None, weight_type='MAP_WEIGHT', 
                     gain='header', texp='header', ncoadd=1, satur=65535., 
                     catformat='ldac', det_area=10, det_thresh=1.5,
                     logfile=None, verbose=True):
    """
    Calls make_fits_cat, but for a generic input fits file
    """

    """ Set output file if user wants default values """
    if outcat=='default':
        outcat = fitsfile.replace('.fits', '.cat')
    if regfile=='default':
        regfile = fitsfile.replace('.fits', '.reg')
    if whtfile=='default':
        whtfile = fitsfile.replace('.fits', '_wht.fits')

    """ Set up for reading information from input file, if requested """
    f = pf.open(fitsfile)
    hdr = f[0].header

    """ Set exposure time """
    if texp == 'header':
        readok = True
        try:
            texp = hdr['exptime']
        except:
            texp = 1.
            readok = False
        if readok:
            print "Exposure time from fits header: %8.1f sec" % texp
        else:
            print "Failed to read EXPTIME header. Setting texp = 1 sec"
    else:
        print "Exposure time set by function call to %8.1f sec" % texp

    """ Set gain """
    if gain == 'header':
        readok = True
        try:
            gain = hdr['gain']
        except:
            gain = 1.
            readok = False
        if readok:
            print "Gain from fits header: %6.3f " % gain
        else:
            print "Failed to read GAIN header. Setting gain = 1.0"
    else:
        print "Gain set by function call to %6.3f" % gain


    f.close()

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur,
                  catformat=catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# ---------------------------------------------------------------------------

def make_fits_cat(fitsfile, outcat='tmp.cat', configfile='sext_astfile.config',
                  gain=1.0, texp=1.0, ncoadd=1, satur=64000., zeropt=None,
                  catformat='ldac', det_area=-1, det_thresh=-1., seeing=0.0, 
                  whtfile=None, weight_type='MAP_WEIGHT', weight_thresh=None,
                  regfile=None, flag_file=None, logfile=None, verbose=True,
                  racol=None, deccol=None, fluxcol=None,
                  fluxerrcol='fluxerr_auto'):
    """
    Runs SExtractor to create the LDAC catalog for the file that the astrometry
    will be applied to.

    Inputs:
        fitsfile    -  input fits file
        configfile -  input SExtractor configuration file.  Note that the
                          gain parameter in this file will be overridden based
                          on the values of the gain, texp, and ncoadd parameters
        gain         -  native gain, in e-/ADU of the instrument
        texp         -  exposure time of the data.  NOTE: keep this at the default
                          value if the units of the data are e- or ADU.  Only change
                          the value from the default if the data are in e-/sec
        ncoadd      -  number of images added together to produce the data in
                          the input fits file.  NOTE: Only change this value if
                          the data were averaged to create the final image.  
                          Keep at the default value, ncoadd=1 if the data were summed.
                          Change to 2xN/3 if the median of N exposures were used.
        satur        -  The chip saturation level, in ADU.  NOTE - this will be
                          multiplied by ncoadd before being passed to SExtractor
        catformat  -  format for the output file.  
                          Possible values: 'ldac', 'ascii'
                          The default, 'ldac', is what is needed for scamp
                          'ascii' produces a catalog that is easier to examine by eye
        det_area    -  Minimum number of pixels for a SExtractor detection, if
                          the value in the configuration file needs to be overridden.
    """

    """ Set up the gain and format of the output catalog (default is LDAC) """
    if catformat.lower() == 'ascii':
        cattype = 'ASCII_HEAD'
    else:
        cattype = 'FITS_LDAC'
    gain_eff = gain * texp * ncoadd

    """ Prepare other parameters """
    satur_eff = satur * ncoadd

    """ Prepare the optional portions of the call to SExtractor """
    sopts = '-GAIN %f -CATALOG_NAME %s -CATALOG_TYPE %s -SATUR_LEVEL %f ' %\
         (gain_eff, outcat, cattype, satur_eff)
    if det_area>0:
        sopts += '-DETECT_MINAREA %d ' % det_area
    if det_thresh>0.:
        sopts += '-DETECT_THRESH %5.2f -ANALYSIS_THRESH %5.2f ' \
             % (det_thresh, det_thresh)
    if zeropt is not None:
        sopts += '-MAG_ZEROPOINT %5.2f ' % zeropt
    if seeing>0.0:
        sopts += '-SEEING_FWHM %6.3f ' % seeing
    if whtfile is None:
        sopts += '-WEIGHT_TYPE NONE '
    else:        
        sopts += '-WEIGHT_TYPE %s -WEIGHT_IMAGE %s ' % (weight_type, whtfile)
    if weight_thresh is not None:
        sopts += '-WEIGHT_THRESH %9.2f ' % weight_thresh
    if flag_file is not None:
        sopts += '-FLAG_IMAGE %s ' % flag_file
    if verbose is False:
        sopts += '-VERBOSE_TYPE QUIET '

    """ Run SExtractor """
    if verbose:
        print ""
        print "Running SExtractor on %s" % fitsfile
        print "Configuration file: %s" % configfile
        print "Override variables:"
        print sopts
        print ""
    if logfile is None:
        try:
            os.system('sex -c %s %s %s' % (configfile, fitsfile, sopts))
        except:
            print ""
            print "ERROR.  Could not run SExtractor on %s" % fitsfile
            print ""
            return
    else:
        try:
            os.system('sex -c %s %s %s > %s' % \
                             (configfile, fitsfile, sopts, logfile))
        except:
            print ""
            print "ERROR.  Could not run SExtractor on %s" % fitsfile
            print ""
            return
    if verbose:
        print ""
        print "Ran SExtractor on %s to produce output catalog %s" % \
             (fitsfile, outcat)

    if regfile is not None:
        if verbose:
            print "Creating ds9 regions file %s from SExtractor catalog." % regfile
        if catformat=='ascii':
            if racol is None:
                racol = 1
                namecol = 0
            if deccol is None:
                deccol = 2
                namecol = 0
            tmpcat = cf.Secat(outcat, catformat, racol=racol, deccol=deccol,
                                            namecol=namecol)
        else:
            tmpcat = cf.Secat(outcat, catformat)
        tmpcat.make_reg_file(regfile, fluxcol, fluxerrcol)
        #make_reg_file(outcat, regfile, catformat)
    if verbose:
        print ""

# -----------------------------------------------------------------------

def make_astrom_cat(fitsfile, configfile='sext_astfile.config', 
                    outcat='tmp.cat', catformat='ldac'):
    """
    Creates an astrometric catalog from a fits file.  This catalog is then
    used as the reference catalog when running scamp.  The catalog will
    be saved in fits LDAC format.
    """

    """ Set up the format of the output catalog (default is LDAC) """
    if catformat.lower() == 'ascii':
        cattype = 'ASCII_HEAD'
    else:
        cattype = 'FITS_LDAC'

    """ Run SExtractor """
    try:
        os.system('sextractor -c %s %s -CATALOG_NAME %s -CATALOG_TYPE %s' % \
                         (configfile, fitsfile, outcat, cattype))
    except:
        print ""
        print "ERROR.  Could not run SExtractor on %s" % fitsfile
        print ""

# -----------------------------------------------------------------------

def make_cat_irac(fitsfile, outcat='tmp.cat', regfile=None,
                  configfile='sext_astfile.config', 
                  gain=29., texp=1., ncoadd=1, satur=64000., catformat='ldac',
                  whtfile=None, weight_type='MAP_WEIGHT', 
                  det_area=15, det_thresh=2.5, logfile=None):
    """
    Calls make_fits_cat, but with gain preset for IRAC channel 2 (4.5 micron).
    Note that the value in the gain parameter is equal to the conversion
    factor given for e-/s to MJy/sr
    """

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur,
                  None, catformat, whtfile=whtfile, weight_type=weight_type,
                  regfile=regfile,
                  det_area=det_area, det_thresh=det_thresh, logfile=logfile)

# -----------------------------------------------------------------------

def make_cat_acs(fitsfile, outcat='default', regfile='default',
                 configfile='sext_astfile.config', 
                 gain=2.0, texp='header', ncoadd=1, satur=65535., det_area=30,
                 det_thresh=2.5, seeing=0.105, obsfilt=None, magsys='ab', 
                 whtfile=None, weight_type='MAP_RMS', 
                 catformat='ldac', logfile=None, verbose=True):
    """
    Calls make_fits_cat, but with parameters preset for HST/ACS WFC
    """

    """ Set output file if user wants default values """
    if outcat=='default':
        outcat = fitsfile.replace('.fits', '.cat')
    if regfile=='default':
        regfile = fitsfile.replace('.fits', '.reg')

    """ Set zeropoint """
    if obsfilt is not None:
        if obsfilt.lower()=='f814w':
            if magsys.lower()=='ab':
                magzp = 25.943
            else:
                magzp = 25.520
    else:
        """ Default zeropoint """
        magzp = 30.

    """ Set exposure time """
    f = pf.open(fitsfile)
    hdr = f[0].header
    if texp == 'header':
        readok = True
        try:
            texp = hdr['exptime']
        except:
            texp = 1.
            readok = False
        if readok:
            print "Exposure time from fits header: %8.1f sec" % texp
        else:
            print "Failed to read EXPTIME header. Setting texp = 1 sec"
    else:
        print "Exposure time set by function call to %8.1f sec" % texp
    f.close()

    """ Call SExtractor """
    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, magzp,
                  catformat, seeing=seeing,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_wfc3(fitsfile, outcat='default', regfile='default',
                  configfile='sext_astfile.config', 
                  whtfile=None, weight_type='MAP_WEIGHT', 
                  gain=2.5, texp='header', ncoadd=1, satur=65535., 
                  catformat='ldac', det_area=30, det_thresh=2.5,
                  logfile=None, verbose=True):
    """
    Calls make_fits_cat, but with gain preset for WFC3
    """

    """ Set output file if user wants default values """
    if outcat=='default':
        outcat = fitsfile.replace('.fits', '.cat')
    if regfile=='default':
        regfile = fitsfile.replace('.fits', '.reg')

    """ Set exposure time """
    f = pf.open(fitsfile)
    hdr = f[0].header
    if texp == 'header':
        readok = True
        try:
            texp = hdr['exptime']
        except:
            texp = 1.
            readok = False
        if readok:
            print "Exposure time from fits header: %8.1f sec" % texp
        else:
            print "Failed to read EXPTIME header. Setting texp = 1 sec"
    else:
        print "Exposure time set by function call to %8.1f sec" % texp
    f.close()

    """ Call SExtractor """
    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur,
                  catformat=catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_wfpc2(fitsfile, outcat='default', regfile='default',
                   configfile='sext_astfile.config', 
                   whtfile=None, weight_type='MAP_WEIGHT', 
                   gain=2.5, texp='header', ncoadd=1, satur=65535., 
                   catformat='ldac', det_area=30, det_thresh=2.5,
                   logfile=None, verbose=True):
    """
    Calls make_fits_cat, but with gain preset for WFPC2
    """

    """ Set output file if user wants default values """
    if outcat=='default':
        outcat = fitsfile.replace('.fits', '.cat')
    if regfile=='default':
        regfile = fitsfile.replace('.fits', '.reg')

    """ Set exposure time """
    f = pf.open(fitsfile)
    hdr = f[0].header
    if texp == 'header':
        readok = True
        try:
            texp = hdr['exptime']
        except:
            texp = 1.
            readok = False
        if readok:
            print "Exposure time from fits header: %8.1f sec" % texp
        else:
            print "Failed to read EXPTIME header. Setting texp = 1 sec"
    else:
        print "Exposure time set by function call to %8.1f sec" % texp
    f.close()

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur,
                  catformat=catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_fors2(fitsfile, outcat='tmp.cat', regfile=None,
                   configfile='sext_astfile.config', 
                   gain=1.25, texp=1., ncoadd=1, satur=64000., catformat='ldac',
                   whtfile=None, weight_type='MAP_WEIGHT', logfile=None):
    """
    Calls make_fits_cat, but with gain preset for FORS2
    """

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, none,
                  catformat, regfile=regfile, 
                  whtfile=whtfile, weight_type=weight_type,
                  logfile=logfile)

# -----------------------------------------------------------------------

def make_cat_isaac(fitsfile, outcat='tmp.cat', regfile=None,
                   configfile='sext_astfile.config', 
                   gain=4.5, texp=1., ncoadd=1, satur=40000., catformat='ldac',
                   whtfile=None, weight_type='MAP_WEIGHT', logfile=None):
    """
    Calls make_fits_cat, but with gain preset for VLT/ISAAC
    """

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, none,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type,
                  logfile=logfile, regfile=regfile)

# -----------------------------------------------------------------------

def make_cat_niri(fitsfile, outcat='tmp.cat', regfile=None,
                  configfile='sext_niri.config', gain=12.3, texp=1., 
                  ncoadd=1, satur=200000., zeropt=None, catformat='ldac',
                  whtfile=None, weight_type='MAP_WEIGHT', flag_file=None,
                  det_thresh=-1, det_area=-1, logfile=None, verbose=True):
    """

    Calls make_fits_cat, but with gain preset for Gemini-North/NIRI

    Note that the read noise (so far not set for this function) is 
    35 e-/pix for the "Medium Background" mode, which is the one typically
    used for broad-band NIR imaging.

    """

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, zeropt,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_wirc(fitsfile, outcat='tmp.cat', regfile=None,
                  configfile='sext_astfile.config', gain=5.467, texp=1., 
                  ncoadd=1, satur=50000., zeropt=None, catformat='ascii',
                  whtfile=None, weight_type='MAP_WEIGHT', flag_file=None,
                  det_thresh=-1, det_area=-1, logfile=None, verbose=True):
    """

    Calls make_fits_cat, but with gain preset for Palomar P200/WIRC

    Note that the read noise (so far not set for this function) is 
    12 e-/pix.

    """

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, zeropt,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_moircs(fitsfile, outcat='tmp.cat', regfile=None,
                    configfile='sext_moircs.config', 
                    ncoadd=1, satur=50000., zeropt=None, catformat='ldac',
                    whtfile=None, weight_type='MAP_WEIGHT', flag_file=None,
                    det_thresh=-1, det_area=-1, logfile=None, verbose=True):
    """

    Calls make_fits_cat, but gets gain from FITS header

    Note that the read noise (so far not set for this function) is 
    ~30 e-/pix.

    """

    """ Get gain and exposure time from header """
    hdr = pf.getheader(fitsfile)
    try:
        gain = hdr['gain']
    except:
        gain = 1.0
    try:
        texp = hdr['exptime']
    except:
        texp = 1.0
    if verbose:
        print ""
        print "File: %s has gain=%6.3f and t_exp = %7.1f" % (fitsfile, gain, texp)


    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, zeropt,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_kait(fitsfile, outcat='tmp.cat', regfile=None,
                  configfile='sext_kait.config', 
                  ncoadd=1, satur=50000., zeropt=None, catformat='ldac',
                  whtfile=None, weight_type='MAP_WEIGHT', flag_file=None,
                  det_thresh=-1, det_area=-1, logfile=None, verbose=True):
    """

    Calls make_fits_cat, but gets gain from FITS header

    Note that the read noise (so far not set for this function) is 
    12 e-/pix.

    """

    """ Get gain and exposure time from header """
    hdr = pf.getheader(fitsfile)
    try:
        gain = hdr['ccdgain']
    except:
        gain = 1.0
    try:
        texp = hdr['exptime']
    except:
        texp = 1.0
    if verbose:
        print ""
        print "File: %s has gain=%6.3f and t_exp = %7.1f" % (fitsfile, gain, texp)


    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, zeropt,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)

# -----------------------------------------------------------------------

def make_cat_hawki(fitsfile, outcat='tmp.cat', regfile=None,
                   configfile='sext_astfile.config', 
                   ncoadd=1, satur=50000., zeropt=None, catformat='ldac',
                   whtfile=None, weight_type='MAP_WEIGHT', 
                   flag_file=None, det_thresh=-1, det_area=-1, 
                   logfile=None, verbose=True):
    """
    Calls make_fits_cat, but gets gain first from the fits file

    Note that readnoise for the HAWK-I chips is 12 e-, or 5 e- if doing
    a non-destructive read of more than 10 samples.
    """

    """ Get gain and exposure time from header """
    hdr = pf.getheader(fitsfile)
    try:
        gain = hdr['gain']
    except:
        gain = 1.0
    try:
        texp = hdr['exptime']
    except:
        texp = 1.0
    if zeropt == 'header':
        try:
            zeropt = hdr['magzpt']
        except:
            print ''
            print 'WARNING: MAGZPT not found in header.'
            zeropt = float(raw_input('Enter zero point to use: '))

    if verbose:
        print ""
        print 'Information from file header: %s' % fitsfile
        print '------------------------------------------------------------------'
        print 'Exposure time: %7.1f' % texp
        print 'Gain:             %6.3f' % gain
        if zeropt == 'header':
            print 'Zero point:     %6.3f' % zeropt

    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, zeropt,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)


# -----------------------------------------------------------------------

def make_cat_suprimecam(fitsfile, outcat='tmp.cat', regfile=None,
                        configfile='sext_scam.config', 
                        ncoadd=1, satur=50000., zeropt=None, catformat='ldac',
                        whtfile=None, weight_type='MAP_WEIGHT', 
                        flag_file=None, det_thresh=-1, det_area=-1, 
                        logfile=None, verbose=True):
    """
    Calls make_fits_cat, but gets gain first from the fits file

    Note that readnoise for the SuprimeCam chips is 10 e-
    """

    """ Get gain and exposure time from header """
    hdr = pf.getheader(fitsfile)
    try:
        gain = hdr['gain']
    except:
        gain = 1.0
    try:
        texp = hdr['exptime']
    except:
        texp = 1.0
    if verbose:
        print ""
        print "File: %s has gain=%6.3f and t_exp = %7.1f" % (fitsfile, gain, texp)


    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur, zeropt,
                  catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)


# -----------------------------------------------------------------------

def run_suprimecam_full(inroot, regroot=None, 
                        configfile='sext_suprimecam.config', 
                        ncoadd=1, satur=50000., zeropt=None, catformat='ldac',
                        whtfile=None, weight_type='MAP_WEIGHT', 
                        flag_file=None, det_thresh=-1, det_area=-1, 
                        logfile=None, verbose=True):
    """
    Calls make_cat_suprimecam for each of the 10 SuprimeCam chips in turn.
    This expects the different chips to be in the form produced by SDFRED2,
     i.e.,  [inroot]_[chipname].fits, where chipname is one of the following:
        chihiro
        clarisse
        fio
        kiki
        nausicaa
        ponyo
        san
        satsuki
        sheeta
        sophie
    """

    chipname = ['chihiro', 'clarisse', 'fio', 'kiki', 'nausicaa', 'ponyo', 'san',
                    'satsuki', 'sheeta', 'sophie']

    for i in chipname:
        fitsfile = '%s_%s.fits' % (inroot, i)
        outcat = '%s_%s.cat' % (inroot, i)
        if regroot is not None:
            regfile = '%s_%s.reg' % (regroot, i)
        make_cat_suprimecam(fitsfile, outcat, regfile=regfile, configfile=configfile,
                                  catformat=catformat)
        

# -----------------------------------------------------------------------

def make_cat_vhs(fitsfile, outcat='tmp.cat', regfile=None,
                 configfile='sext_vhs.config', 
                 ncoadd=1, satur=50000., catformat='ldac',
                 whtfile=None, weight_type='MAP_RMS', 
                 flag_file=None, det_thresh=-1, det_area=-1, 
                 logfile=None, verbose=True):
    """
    Calls make_fits_cat, but gets some relevant info first from the fits file
    """

    """ Get information from header """
    hdr = pf.getheader(fitsfile)
    """
    For gain use median value from VISTA/VHS web site
     http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/vista-gain
    """
    try:
        gain = hdr['gain']
    except:
        gain = 4.19
    try:
        texp = hdr['texptime'] # NB: different from most
    except:
        texp = 1.0
    try:
        zeropt = hdr['photzp']
    except:
        print ''
        print 'Warning: PHOTZP not found in header. '
        zeropt = float(raw_input('Enter zero point to use: '))
    if verbose:
        print ""
        print 'Information from file header'
        print '----------------------------'
        print 'Exposure time: %7.1f' % texp
        print 'Zero point:     %6.3f' % zeropt


    make_fits_cat(fitsfile, outcat, configfile, gain, texp, ncoadd, satur,
                  zeropt, catformat,
                  whtfile=whtfile, weight_type=weight_type, 
                  det_thresh=det_thresh, det_area=det_area, flag_file=flag_file,
                  logfile=logfile, regfile=regfile, verbose=verbose)


# -----------------------------------------------------------------------

def remove_nans(infile, outfile, replaceval=999999., hext=0):
    """
    SExtractor does not extract any sources if the fits file contains
    nan values.  Therefore, change the values of any pixel containing
    a nan to something that exceeds the saturation value in the 
    SExtractor config file.

    Inputs:
        infile      - input file that may contain nan values
        outfile     - output file with nans replaced
        replaceval - value to assign to pixels that have nan values
        hext         - HDU containing data
    """

    import imfuncs as im

    """ Get data from input file """
    hdulist = im.open_fits(infile)
    if hdulist is None:
        print ""
        print "ERROR.  Could not open input file %s" % infile
        print ""
        return
    hdr  = hdulist[hext].header.copy()
    data = hdulist[hext].data.copy()

    """ Check for nans and replace any """
    mask = n.isnan(data)
    if mask.sum()>0:
        data[mask] = replaceval

    """ Write output """
    pf.PrimaryHDU(data, hdr).writeto(outfile, clobber=True)
    del data, hdr, hdulist

# -----------------------------------------------------------------------

def run_scamp(incat, asttype='file', astcat=None, 
                  configfile='scamp_reffile.config'):
    """
    Runs scamp on the LDAC catalog produced by running make_fits_cat on the
     fits file for which the astrometric solution is to be calculated.
    The astrometric reference can either be one of the standard catalogs
     (2MASS, USNO-B, etc.) or a catalog produced from a local fits file
     (probably by running make_astrom_cat).

    Inputs:
        incat        - input catalog, produced by make_fits_cat from the fits file 
                         for which astrometric solution is to be calculated
        asttype     - astrometric reference.  For now this can be '2MASS',
                         'USNO-B1', 'SDSS-R5', or 'FILE'
        astcat      - if asttype is FILE, then this parameter identifies the
                         reference astrometric catalog (produced by make_astrom_cat)
        configfile - scamp configuration file

    """

    """ Check the astrometric reference """
    print asttype.upper()
    if asttype.upper()=='2MASS' or asttype.upper()=='USNO-B1' or \
           asttype.upper=='SDSS-R5':
        refcat = asttype.upper()
    elif asttype.upper()=='FILE':
        print "found file"
        refcat = asttype.upper()
        if astcat is None:
            print ""
            print "ERROR. Astrometric reference is FILE but no filename given"
            print "         in the astcat parameter"
            return
    else:
        print ""
        print "ERROR. Astrometric reference type (%s) is not recognized" % asttype
        print "For now, only valid values are 2MASS, USNO-B1, SDSS-R5, or FILE"
        print ""
        return
    if refcat == 'FILE':
        try:
            os.system('scamp %s -c %s -ASTREF_CATALOG FILE -ASTREFCAT_NAME %s' % \
                         (incat, configfile, astcat))
        except:
            print ""
            print "ERROR.  Could not run scamp on %s" % incat
            print ""
    else:
        try:
            os.system('scamp %s -c %s -ASTREF_CATALOG %s' % \
                         (incat, configfile, refcat))
        except:
            print ""
            print "ERROR.  Could not run scamp on %s" % incat
            print ""

# ------------------------------------------------------------------------------

def import_ascii_header(fitsfile, headfile, hext=0):
    """
    If it completes successfully, scamp will produce an ascii representation
    of the WCS header cards.  This function will read the ascii file and
    update the relevant headers in the fits file with the information
    derived from the ascii file.

    Inputs:
        fitsfile  - fits file for which the header information will be updated
        headfile  - ascii header file that has been produced by scamp
        hext        - designates which HDU in the input fits file to update
    """

    """ Read the input fits file """
    import imfuncs as im
    hdulist = im.open_fits(fitsfile, 'update')
    if hdulist is None:
        return
    hdr = hdulist[hext].header

    """ Read the ascii file and update the relevant header cards """
    print ""
    print "Updating/adding WCS header cards for %s" % fitsfile
    print "--------------------------------------------------------------"
    for line in open(headfile, 'r'):
        foo = line.strip().split('=')
        field = foo[0][0:9]
        wcskeys = ['RADEC', 'CTYPE', 'CUNIT']
        for key in wcskeys:
            if field.upper()[0:5] == key:
                tmp = foo[1].split('/')
                val = str(tmp[0].strip()).strip('\'')
                print "%-8s %s" % (field, val)
                if len(tmp) == 2:
                    hdr.update(field, val, tmp[1])
                else:
                    hdr.update(field, val)
        wcskeys = ['CRVAL', 'CRPIX', 'EQUIN']
        for key in wcskeys:
            if field.upper()[0:5] == key:
                tmp = foo[1].split('/')
                val = float(tmp[0])
                print "%-8s %f" % (field, val)
                if len(tmp) == 2:
                    hdr.update(field, val, tmp[1])
                else:
                    hdr.update(field, val)
        wcskeys = ['CD1', 'CD2', 'PV1', 'PV2']
        for key in wcskeys:
            if field.upper()[0:3] == key:
                tmp = foo[1].split('/')
                val = float(tmp[0])
                print "%-8s %f" % (field, val)
                if len(tmp) == 2:
                    hdr.update(field, val, tmp[1])
                else:
                    hdr.update(field, val)

    """ Write out updated fits file and clean up """
    hdulist.flush()
    del hdr, hdulist

# ------------------------------------------------------------------------------

def run_swarp(infile, outfile, configfile='swarp.config', pixscale=None, 
              whttype=None, wht_suffix=None, combtype=None, backsub=None,
              centertype=None, center=None, imsize=None, proj_type=None):
    """
    Runs swarp on the input image (infile) to produce a resampled, and possibly
    coadded output file (outfile).
    ** NOTE: For coaddition, the input "file" should actually be a file 
        containing a list of input files (one per line) and this input filename
        should be preceded by the @ character.  Thus, to coadd three file, one
        might have a file called 'filelist.txt' that contains the following three
        lines:
          file1.fits
          file2.fits
          file3.fits
      In this case, infile would be '@filelist.txt'

    Inputs:
      infile      - input fits file with new WCS in header or file with file list
      outfile     - swarped output file
      configfile - swarp configuration file

     *** The following parameters should only be set if the user wants to
     *** override the values in the configuration file
      pixscale    - pixel scale in arcsec/pix
      whttype     - type of weight file(s) that are associated with the input
                        fits file(s).  Possible values are NONE, MAP_WEIGHT,
                        BACKGROUND, MAP_RMS, or MAP_VARIANCE (see swarp manual)
      wht_suffix - suffix for weight file(s) associated with input file(s), e.g.,
                        '_wht.fits'
      combtype    - method used to coadd input images, if the input file is a
                        list of files.  Possible values are MEDIAN, AVERAGE, MIN, 
                        MAX, WEIGHTED, or CHI2 (sww swarp manual)
      backsub     - Subtract background? Possible values are Y or N
      centertype - Method for determining center.  Possible values are
                        MANUAL, ALL or MOST
      center      - center for projected image, if centertype is MANUAL (see 
                        swarp manual)
    """

    """
    Start by setting parameters to be passed to swarp
    """
    whtout = outfile.replace('.fits', '_wht.fits')
    swopts = '-IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s ' % (outfile, whtout)
    if pixscale is not None:
        swopts += '-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %7.4f ' % pixscale
    if whttype is not None:
        swopts += '-WEIGHT_TYPE %s ' % whttype
    if wht_suffix is not None:
        swopts += '-WEIGHT_SUFFIX %s ' % wht_suffix
    if combtype is not None:
        swopts += '-COMBINE_TYPE %s ' % combtype
    if backsub is not None:
        swopts += '-SUBTRACT_BACK %s ' % backsub
    if centertype is not None:
        swopts += '-CENTER_TYPE %s ' % centertype
    if center is not None:
        swopts += '-CENTER %s ' % center
    if imsize is not None:
        swopts += '-IMAGE_SIZE %d ' % imsize

    """ Run swarp """
    print ""
    print "Running swarp on     %s" % infile
    print "Configuration file: %s" % configfile
    print "Output file:          %s" % outfile
    print "Override variables:"
    print swopts
    print ""
    try:
        os.system('swarp %s -c %s %s' % (infile, configfile, swopts))
    except:
        print ""
        print "ERROR.  Could not run swarp on %s" % infile
        print ""
        return
    print ""
    print "Finished running swarp %s -c %s %s" % (infile, configfile, swopts)

# ------------------------------------------------------------------------------

def do_photom(infiles, photocat, magcol, photzp=30., incats='default',    
              fixzp=None, mautocol=19, mapermaxcol=49, returnzp=False, 
              racol=1, deccol=2, xcol_fits=8, ycol_fits=9, max_offset=40.,
              phottype='2MASS', verbose=True, summary=True, doplot=False):
    """

    NOT YET IMPLEMENTED

    Use a photometric catalog to set the proper zero-point for the input fits 
    file(s).  This function will run SExtractor on the input fits 
    file(s) and get the instrumental magnitude of any of the photometric objects 
    in the field of view of the image.  By comparison to the true magnitudes, the
    median zero-point for the image will be determined.
    """

    """ Set up for running the task """
    if type(infiles) is str:
        if verbose:
            print ""
            print "Single input file"
        tmplist = [infiles,]
        tmpcat = [incats,]
    elif type(infiles) is list:
        if verbose:
            print ""
            print "Input file list"
        tmplist = infiles
        if type(incats) is str:
            tmpcat = []
            for i in range(len(infiles)):
                tmpcat.append(incats)
        else:
            tmpcat  = incats
    else:
        print ""
        print "Warning.  Input files need to be either a list of files "
        print " (python type==list) or a single input file name."
        print ""
        return

    """ Get the data from the photometric catalog """
    photdat = n.loadtxt(photocat)

    """ Loop over the input file(s) """
    zpmedind  = n.zeros(len(tmplist))
    count = 0
    for i in range(len(tmplist)):
        """ Open up the input file and extract information """
        f = tmplist[i]
        hdu = pf.open(f, mode='update')
        hdr = hdu[0].header
        if tmpcat[i] == 'default':
            fitscat = f.replace('.fits', '_phot.cat')
        else:
            fitscat = tmpcat[i]

        """ Find matches between the image and the photometric catalog """
        secat = cf.Secat(fitscat, verbose=verbose)
        secat.match_fits_to_ast(f, photocat, max_offset=max_offset, racol=racol,
                                deccol=deccol, xcol=xcol_fits, ycol=ycol_fits,
                                doplot=doplot, verbose=True)
        tmpphot = photdat[secat.astmask][secat.goodmask]
        tmpdat  = secat.data[secat.matchind]
        dauto = tmpdat['f%d' %mautocol] - tmpphot[:, magcol]
        daper = tmpdat['f%d' %mapermaxcol] - tmpphot[:, magcol]
        dmmed_auto = n.median(dauto)
        dmmed_aper = n.median(daper)
        if secat.nmatch>2:
            dmsig_auto = dauto.std()
            dmsig_aper = daper.std()
        else:
            dmsig_auto = -99.
            dmsig_aper = -99.
        if verbose:
            print ''
            print 'Matched %d objects to %s stars' % (secat.nmatch, phottype)
            print '--------------------------------------'
            print '  %s ID m_%s mauto maper  dmauto dmaper' % (phottype, phottype)
            print ' --------- ------- ----- -----  ------ ------'
            for i in range(secat.nmatch):
                print ' %9.0f  %5.2f  %5.2f %5.2f  %6.2f %6.2f' % \
                    (tmpphot[i,0], tmpphot[i, magcol], tmpdat[i]['f%d'%mautocol],
                     tmpdat[i]['f%d'%mapermaxcol], dauto[i], daper[i])
            print ''
            print 'dm_auto = %6.2f +/- %6.2f.    dm_aper = %6.2f +/- %6.2f' % \
                 (dmmed_auto, dmsig_auto, dmmed_aper, dmsig_aper)
            print 'zp_auto = %6.2f +/- %6.2f.    zp_aper = %6.2f +/- %6.2f' % \
                 (photzp-dmmed_auto, dmsig_auto, photzp-dmmed_aper, dmsig_aper)

        """ Store information on the photometry in the fits header """
        hdr.update('nmatchph', secat.nmatch, 'Number of %s stars in image' %
                      phottype)
        hdr.update('magzp', photzp-dmmed_aper,
                      'Magnitude zeropoint relative to %s' % phottype)
        hdr.update('magzpsig', dmsig_aper, 'Magnitude zeropoint uncertainty')
        hdu.flush()
        zpmedind[count] = photzp - dmmed_aper
        count += 1

    """ 
    Further modify the fits header to contain the FLXSCALE keyword that is
    used by swarp.
    """
    if summary:
        print ''
        print 'File                 N*     ZP  sig_ZP FLXSCALE'
        print '---------------  --  ----- ------ --------'
    for i in range(len(tmplist)):
        hdu = pf.open(tmplist[i], mode='update')
        f = tmplist[i].replace('.fits', '')
        hdr = hdu[0].header
        if fixzp is not None:
            zpfid = fixzp
        else:
            zpfid = n.median(zpmedind)
        fluxscale = 10.**((hdr['magzp'] - zpfid)/-2.5)
        hdr.update('flxscale', fluxscale, 'Scale factor needed for swarp coaddition')
        if summary:
            print '%-15s  %2d  %5.2f %6.2f  %6.4f' % \
                 (f, hdr['nmatchph'], hdr['magzp'], hdr['magzpsig'], hdr['flxscale'])
        hdu.flush()
        del hdr
    if summary:
        print ''
        print 'Median zero-point for these files:  %5.2f' % n.median(zpmedind)
        print 'FLXSCALE values set for zero-point: %5.2f' % zpfid
        print ' (different values => FLXSCALE set to a fixed value set by'
        print '  the passed fixzp parameter)'
        print ''
    
    if returnzp:
        return n.median(zpmedind)
