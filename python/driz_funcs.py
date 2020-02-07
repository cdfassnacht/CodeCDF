import os
import glob
import shutil
import numpy as np
from math import sqrt

from astropy import units as u
from astropy.io import fits as pf
from astropy.io import ascii
from astropy.table import Table

from drizzlepac import tweakreg
from drizzlepac import tweakback
from drizzlepac import astrodrizzle as astd

from specim.imfuncs import image as imf

# ---------------------------------------------------------------------------

class DrizSet():
    """

    A class set up to make some astrodrizzle-related processing easier

    """

    def __init__(self, archdir, archstr, inst, skystat=True, verbose=True):
        """

        Sets up the lists of files, and reports some statistics if desired

        Required inputs:
          archdir - directory containing the *flt.fits or *flc.fits that
                    that were downloaded from the HST archive
          archstr - string (with wildcard symbols) that defines the 
                    the files in the archive directory.
                    Examples: '*.fits', 'i*flc.fits', 'icv*flc.fits'
          inst    - instrument used for the observation.  Current options are
                     acs, wfc3uv, wfc3ir

        Optional inputs:
          skystat - Set to True (the default) to get some statistics for
                    the input files
        """

        self.inst = inst
        
        self.archdir = archdir
        self.archfiles = glob.glob('%s/%s' % (archdir, archstr))
        self.archfiles.sort()
        self.nfiles = len(self.archfiles)
        self.infiles = []
        for f in self.archfiles:
            self.infiles.append(os.path.basename(f))

        if skystat:
            self.sky_stats(verbose)

    # -----------------------------------------------------------------------

    def sky_stats(self, verbose=True):
        """
        Reports several image statistics for the input images, including
         sky level  - from clipped mean
         sky rms    - from clipped rms, in e-
        """

        """
        Set edge region to exclude for statistics, via the width in pixels
        of the exlusion region along each edge of the detector
        """
        wexcl = 30
    
        """ Set up table to hold the information """
        fname = np.zeros(self.nfiles, dtype='S15')
        texp = np.zeros(self.nfiles)
        mean1 = np.zeros(self.nfiles)
        rms1 = np.zeros(self.nfiles)
        rmsps1 = np.zeros(self.nfiles)
        mean2 = np.zeros(self.nfiles)
        rms2 = np.zeros(self.nfiles)
        rmsps2 = np.zeros(self.nfiles)
        tab = Table([fname, texp, mean1, rms1, mean2, rms2],
                    names=('fname', 'texp', 'mean1', 'rms1', 'mean2', 'rms2'))

        """ Set column units and formats """
        tab['fname'].format = '%-15s'
        tab['texp'].format = '%6.1f'
        tab['mean1'].format = '%5.2f'
        tab['rms1'].format = '%4.2f'
        tab['mean2'].format = '%5.2f'
        tab['rms2'].format = '%4.2f'
        tab['texp'].unit = u.s
        tab['mean1'].unit = u.electron
        tab['rms1'].unit = u.electron
        tab['mean2'].unit = u.electron
        tab['rms2'].unit = u.electron
    
        """ Loop through the files"""
        if verbose:
            print('')
            print('Estimating sky statistics')
            print('-------------------------')
        for i, f  in enumerate(self.archfiles):
            tab['fname'][i] = (os.path.splitext(self.infiles[i]))[0]
            if verbose:
                print(tab['fname'][i])
            hdr = pf.getheader(f)
            tab['texp'][i] = hdr['exptime']
            sci1 = imf.Image(f, hext=1, verbose=False, wcsverb=False)
            sci2 = imf.Image(f, hext=4, verbose=False, wcsverb=False)
            in1 = sci1['input']
            in2 = sci2['input']
            ymax, xmax = in1.data.shape
            statsec1 = [wexcl, wexcl, xmax-wexcl, ymax-wexcl]
            ymax, xmax = in2.data.shape
            statsec2 = [wexcl, wexcl, xmax-wexcl, ymax-wexcl]
            in1.sigma_clip(statsec=statsec1)
            in2.sigma_clip(statsec=statsec2)
            tab['mean1'][i] = in1.mean_clip
            tab['rms1'][i] = in1.rms_clip
            tab['mean2'][i] = in2.mean_clip
            tab['rms2'][i] = in2.rms_clip
            del hdr, sci1, sci2
        
        """ Print out the results """
        if verbose:
            print('')
            print(tab)
            print('')
            varsum = (tab['rms1']**2).sum()
            ttot = tab['texp'].sum()
            print('Predicted rms in e-/sec: %6.4f' % (sqrt(varsum)/ttot))

        self.instats = tab
        
    # -----------------------------------------------------------------------

    def driz_prep(self, drizdir):
        """

        Prepares the directory for drizzling

        """

        """ Get rid of old version of drizzle working directory """
        os.system('rm -rf %s' % drizdir)
        os.system('mkdir %s' % drizdir)

        """ Copy the *flc or *flt files into the drizzle directory """
        for f in self.archfiles:
            os.system('cp -v %s %s/.' % (f, drizdir))

    # -----------------------------------------------------------------------

    def set_defaults(self):
        """

        Set the default astrodrizzle parameters, including some that are
        instrument-dependent.

        Note that astrodrizzle itself has default values for all of these
         parameters.  The values below may match the astrodrizzle defaults or
         they may overwrite them.
        In turn, all of the parameters that are set here can be overwritten
         via the keyword arguments (**kwargs) in the run method below.

        The defaults are set in different blocks corresponding to the different
         astrodrizzle steps, except for a few parameters that depend on the
         instrument or number of exposures.  Those defaults are set after the
         generic defaults.

        """

        """ Set the general / initialization defaults """
        defpars = {}
        defpars['updatewcs'] = False
        defpars['clean'] = True

        """ Static mask defaults """
        defpars['static'] = True

        """ Sky subtraction defaults """
        defpars['skysub'] = True

        """ Initial drizzling defaults """
        defpars['driz_separate'] = True
        defpars['driz_sep_kernel'] = 'turbo'
        defpars['driz_sep_scale'] = None
        defpars['driz_sep_pixfrac'] = 1.
        defpars['driz_sep_bits'] = '16, 64'

        """
        Defaults for creating a median image
        NOTE: the combine_type, combine_nhigh, and combine_nlow parameters are
        set below depending on the number of exposures
        """
        defpars['median'] = True

        """
        Set default values for parameters that are sensitive to the number of
        input images
        """
        if self.nfiles > 4:
            defpars['combine_type'] = 'median'
            if self.nfiles % 2 == 0:
                defpars['combine_nhigh'] = 1
            else:
                defpars['combine_nhigh'] = 0
            if self.nfiles > 7:
                defpars['combine_nhigh'] += 2
            if self.nfiles > 11:
                defpars['combine_nhigh'] += 2
        else:
            defpars['combine_type'] = 'minmed'
            defpars['combine_nhigh'] = 0

        """ Blot defaults """
        defpars['blot'] = True

        """ Cosmic ray mask defaults """
        defpars['driz_cr'] = True
        defpars['driz_cr_corr'] = True

        """
        Defaults for final drizzle
        NOTE: the default final scales and pixfrac values are
         instrument-dependent and are set below
        """
        defpars['driz_combine'] = True
        defpars['build'] = False
        defpars['final_wht_type'] = 'EXP'
        defpars['final_kernel'] = 'lanczos3'
        defpars['final_wt_scl'] = 'exptime'
        defpars['final_bits'] = '16, 64'
        defpars['final_wcs'] = True
        defpars['final_units'] = 'cps'
        defpars['final_rot'] = 0

        """ Set the instrument-specific defaults """
        if self.inst.lower() == 'acs':
            defpars['final_scale'] = 0.05
            defpars['final_pixfrac'] = 1.
        elif self.inst.lower() == 'wfc3uv' or inst.lower() == 'wfc3uvis':
            defpars['final_scale'] = 0.04
            defpars['final_pixfrac'] = 1.
        elif self.inst.lower() == 'wfc3ir':
            defpars['final_scale'] = 0.06
            defpars['final_pixfrac'] = 0.7
        else:
            raise ValueError('Instrument parameters not yet implement for %s'
                             % self.inst)

        """ Return the default settings """
        return defpars

    # -----------------------------------------------------------------------

    def call_astrodriz(self, outroot, drizdir, **kwargs):
        """

        Does the actual call to astrodrizzle in the desired directory

        """

        os.chdir(drizdir)
        astd.AstroDrizzle(self.infiles, output=outroot, **kwargs)
        if os.path.isfile('%s_drc_sci.fits' % outroot):
            os.system('cp %s_drc_sci.fits ../%s.fits' % (outroot, outroot))
        for ext in ['ctx', 'wht']:
            if os.path.isfile('%s_drc_%s.fits' % (outroot, ext)):
                os.system('cp %s_drc_%s.fits ../%s_%s.fits' %
                          (outroot, ext, outroot, ext))
        os.chdir('..')

    # -----------------------------------------------------------------------

    def run_adriz(self, outroot, drizdir='finalDrz', debug=False, **kwargs):
        """

        Prepares for the astrodrizzle call and then calls astrodrizzle

        """

        """ Prepare for the call """
        self.driz_prep(drizdir)

        """ Set up the default parameters """
        adpars = self.set_defaults()
        if debug:
            print(adpars)

        """ Now overwrite / append any additional parameters in kwargs """
        if len(kwargs) > 0:
            for k in kwargs.keys():
                if debug:
                    print('kwargs: Setting %s to %s' % (k, str(kwargs[k])))
                adpars[k] = kwargs[k]

        if debug:
            print(adpars)

        """ Call the astrodrizzle routine """
        self.call_astrodriz(outroot, drizdir, **adpars)
    
    # -----------------------------------------------------------------------

    def put_blt_in_orig(self, blotdir, workdir='aligned_files', verbose=True):
        """

        Creates temporary files that are copies of the WCS information of the
         original files from the archive (either *flc.fits or *flt.fits) but
         have their data replaced by the blotted version that are created
         after running astrodrizzle on the original files.
        The motivation for doing this is that running tweakreg is (perhaps)
         easier with images that have been cleaned of all of their cosmic rays.
        This process consists of the following steps:
         1. Copy the original *flc or *flt files into to the output directory
         2. Copy the data from the blot files into the corresponding temporary
            files 

        Required inputs:
          blotdir  - directory containing the blotted files

        Optional inputs:
          workdir  - directory containing the output images.  Default is a
                      directory called 'aligned_files'
          verbose  - set to True (the default) to print out the filenames

        """
    
        for f,t in zip(self.archfiles, self.infiles):
            """ Copy the original files into the working directory """
            os.system('cp -v %s %s' % (f, workdir))
            outfile = os.path.join(workdir, t)
        
            """ Define the files that contain the blotted data """
            base = t.split('_')[0]
            bltfile1 = os.path.join(blotdir, '%s_sci1_blt.fits' % base)
            bltfile2 = bltfile1.replace('sci1', 'sci2')
            if verbose:
                print(t, bltfile1, bltfile2)

            """ Copy the blotted data into the new output files """
            flcblt = pf.open(outfile, mode='update')
            flcblt['sci', 1].data = pf.getdata(bltfile1)
            flcblt['sci', 2].data = pf.getdata(bltfile2)
            flcblt.flush()

    # -----------------------------------------------------------------------

    def run_tweakreg(self, workdir, conv_width, detthresh, searchrad, fitgeom,
                     outshifts, wcsname, updatehdr=False, excludefile=None,
                     minobj=15, refcat=None, **kwargs):
        """

        Calls tweakreg, after setting up some default values that seem to work
         well with this data set

        Required inputs
          workdir     - the working directory, which contains the working
                        version of the files (the ones that contain the
                        blotted data)
          conv_width  - Width of a PSF-like object, used to choose objects that
                         will be used to determine the alignment.  For HST
                         the "default" value is 3.5 (pixels), but this can be
                         larger (e.g., 6) if the field doesn't have a lot of 
                         stars.  The larger value would allow, e.g., compact
                         galaxies to be used as well.
          detthresh   - Detection threshold in terms of sigma above mean.  The
                         tweakreg default is 4.0, but this may find too many
                         objects, many of which may be spurious.  Choose larger
                         values for fewer objects.
          searchrad   - Search radius, in arcsec.  The WCS of the input images
                         are used to predict where the object would be in the
                         reference image.  HST observations obtained at
                         different epoch may have WCS solutions that differ by
                         up to an arcsecond (possibly more in really bad
                         cases) due to different guide stars.
          fitgeom     - Fitting geometry.  Choices are 'rscale' or 'general'
          outshifts   - Name of the output file containing the shifts
          wcsname     - Name for the WCS solution that will be put into the 
                        fits header once the updatehdr parameter is set to True.
                        The tweakreg default may be 'TWEAK' but that is not too
                        informative.

        Optional inputs
          updatehdr   - Should the WCS information in the fits file header be
                         updated?  The default (False) says no.  Only set this
                         to True once you are satisfied that you have found a
                         good WCS solution (by looking at the output shifts
                         table and residual and vector plots)
          excludefile - A file that lists, for each input image, the files that
                         indicate the good and bad regions of the detectors.
                         The object detection in tweakreg will ignore objects in
                         the bad regions.  The default value of None means that
                         all regions of all detectors will be used for object
                         detection.
          minobj      - Minimum number of matched objects to use for a valid
                         solution.  Default value is 15
          refcat      - External reference catalog (e.g., from SDSS).  If this
                         is given (the default value is None), then the images
                         will be compared to the catalog rather than to the
                         reference image
          **kwargs    - Additional tweakreg parameters (see tweakreg help)

        """

        """ Set up input files for tweakreg """
        infiles = []
        for f in self.infiles:
            infiles.append(os.path.join(workdir, f))
        
        """ Run tweakreg """
        tweakreg.TweakReg(infiles, searchrad=searchrad, minobj=minobj,
                          imagefindcfg={'threshold': detthresh,
                                        'conv_width': conv_width},
                          refcat=refcat, interactive=False, shiftfile=True, 
                          outshifts=outshifts, wcsname=wcsname, reusename=True,
                          fitgeometry=fitgeom, exclusions=excludefile,
                          updatehdr=updatehdr, **kwargs)

        """ Check the output rms values """
        shift_tab = ascii.read(outshifts,format='no_header',
                               names=['file', 'dx', 'dy', 'rot', 'scale',
                                      'xrms', 'yrms'])

        """ Add columns for shifts in terms of rms """
        shift_tab['fdx'] = shift_tab['dx'] / shift_tab['xrms']
        shift_tab['fdy'] = shift_tab['dy'] / shift_tab['yrms']

        """ Print out the information """
        formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f', '.4f', '.4f']
        for i, col in enumerate(shift_tab.colnames[1:]):
            shift_tab[col].format = formats[i]
        print(shift_tab)

    # -----------------------------------------------------------------------

    def put_orig_in_blt(self, workdir, debug=False):
        """

        Copies the data from the original *flc / *flt files to the working
         versions in the working directory, which were copies of the original
         *flc.fits / *flt.fits files but with the data replaced by the blotted
         data.

        This step is necessarily run after running the put_blt_in_orig
         function and also after running run_tweakreg

        The motivation is that the working versions, with the clean data from
         the blotting, are used to determine what the updates to the WCS should
         be, via tweakreg.  Now that the WCS has been corrected (in a separate
         tweakreg step), copy the original data back into these files so that
         they can be used to run astrodrizzle again.

        Required inputs:
         workdir - the working directory, which contains the working version
                    of the files (the ones that contain the blotted data)

        Optional inputs:
         debug     - Set to True to get additional info. Default value is False

        """
        for f, t in zip(self.archfiles, self.infiles):
            flcorig = pf.open(f)
            workfile = os.path.join(workdir, t)
            flcblt = pf.open(workfile, mode='update')
            print('Copying data from %s to %s' % (f, t))
            if debug:
                print(flcorig['sci', 1].data.shape)
                print(flcblt['sci', 1].data.shape)
            flcblt['sci', 1].data = flcorig['sci', 1].data.copy()
            flcblt['sci', 2].data = flcorig['sci', 2].data.copy()
            flcblt.flush()

