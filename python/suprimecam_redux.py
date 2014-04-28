"""
A set of functions that can be used for the latter part of the reduction
of Subaru SuprimeCam data.  

Note that the current approach is to use Subaru's SDFRED pipeline for the
following steps:
  File renaming:         namechange.csh
  Overscan subtraction:  overscansub.csh
  Flat-field creation:   mask_mkflat_HA.csh
  Flat-field correction: ffield.csh
  Distortion correction: distcorr.csh
  AG masking:            mask_AGX.csh

The functions in this file take over after the AG masking step.
"""

import pyfits as pf
import numpy as n

#---------------------------------------------------------------------------

def rename_before_swarp(indir='../Calib'):
    """

    Creates symbolic links from to post AG masking step in order to have
      more compact filenames.
    Input files are expected to have the following form:
        [indir]/AgfTo_RH[obsdate]objectnnn_[chipname].fits
    They will be link with a link name of [outdir]/objectnnn_[chipname].fits

    Inputs:
      indir  -  Location of post AG masking files (AgfTo_RH*fits)

    """
    import glob
    import os

    """ Get the input file list """
    infiles = glob.glob('%s/AgfTo_RH*fits' % indir)

    """ Rename the files """
    for f in infiles:
        objchip = f.split('object')[1]
        outfile = 'object%s' % objchip
        os.system('ln -s %s %s' % (f,outfile))
        print 'Linked %s to %s' %(f,outfile)

#---------------------------------------------------------------------------

def make_wht_for_swarp(infiles, mingood=-100, outext='_wht'):
    """
    Creates a weight file for each input file, in preparation for running
    swarp the first time.  The SuprimeCam pipeline marks bad pixels with
    a value of -2^15 = -32768.  However, there are other bad pixels with
    slightly different values, perhaps because of the flat-fielding.
    Therefore, create the weight file using the following algorithm:
    
       wht=1  for all pixels with values => mingood
       wht=0  for all pixels with values < mingood

    Inputs:
       infiles  - list of input files, perhaps created with a glob.glob call
       mingood  - minimum good value.  Default = -100
       outext   - extension used for output weight file name. In other
                    words, for an input file of [root].fits the output file
                    name will be [root][outext].fits
                  SExtractor, swarp, etc., have a default of:   '.weight'
                  However, the default for this function is:    '_wht'

    Outputs:
       Each input file called [root].fits will produce an output weight file
        called [root][outext].fits
    """

    """ Make sure that the input is either a list or a single file """

    if type(infiles) is str:
        print ""
        print "Single input file"
        tmplist = [infiles,]
    elif type(infiles) is list:
        print ""
        print "Input file list with %d members" % (len(infiles))
        tmplist = infiles
    else:
        print ""
        print "Warning.  Input frames need to be either a list of files "
        print " (python type==list) or a single input file name."
        print ""
        return

    """ Run through the list """

    print ''
    print 'Making weight files'
    print '-------------------'
    for f in tmplist:
        """ Open input file and get the object name """
        print 'Input file:  %s' % f
        hdu = pf.open(f)
        data = hdu[0].data
        hdr = hdu[0].header
        try:
            objname = hdr['object']
        except:
            objname = 'Unknown object'

        """ Create the weight data """
        whtdat = n.ones((data.shape[0],data.shape[1]))
        whtdat[data<mingood] = 0

        """ Write the output file """
        ohdu = pf.PrimaryHDU(whtdat)
        ohdr = ohdu.header
        oname = 'Weight file for %s' % objname
        ohdr.update('object',oname,'Object name')
        ofile = f.replace('.fits','%s.fits') % outext
        ohdu.writeto(ofile)
        print 'Output file: %s' % ofile
        print ''

        """ Clean up """
        hdu.close()
        del ohdu

#---------------------------------------------------------------------------

def robust_sigma(data, refzero=False):
    """
    Calculate a robust estimate of the dispersion of a distribution.
    For an uncontaminated distribution, this estimate is identical to
    the standard deviation.

    This code is ported from robust_sigma.pro in IDL/astrolib

    Inputs:

    Output:
       rsig - robust sigma estimator.  Return -1 if failure.
    """

    """ Set a tolerance """
    eps = 1.e-20

    """ Set central point for cfomputing the dispersion """
    if refzero:
        dat0 = 0.0
    else:
        dat0 = n.median(data)
    datdiff = (data - dat0).flatten()

    """ Find absolute deviation about the median """
    mad = n.median(n.absolute(datdiff))/0.6745

    """ Try the mean absolute deviation if the mad is zero """
    if mad<eps:
        mad = (n.absolute(datdiff)).mean()/0.8
    if mad<eps:
        return 0.0

    """ Do the biweighted value """
    u = datdiff / (6. * mad)
    uu = u*u
    q = uu<1.
    if q.sum()<3:
        print ''
        print 'robust_sigma: input distribution is just too weird.'
        print 'returning value of -1.'
        return -1.
    ntot = data[n.isfinite(data)].sum()
    num = ((data[q] - dat0)**2 * (1.-uu[q])**4).sum()
    denom = ((1.-uu[q]) * (1. - 5.*uu[q])).sum()
    rvar = ntot * num / (denom * (denom - 1.))

    if rvar>0.:
        return n.sqrt(rvar)
    else:
        return 0.

#---------------------------------------------------------------------------

def make_wht_for_final(infiles, medfits, sig=4.):
    """
    Creates a weight file for each input image to a final swarp call.
    This weight file will assign zero weight to pixels that differ by
    more than sig*rms from the median-stacked image that was created
    in an earlier step.
    Therefore, this function is a lot like the "blot" step in multidrizzle.

    Inputs:
       infiles - list of individual input fits files that will be compared to the
                 median stack.  These will probably be called *resamp.fits
       medfits - median-stacked fits file
    """

    """ Make sure that the input is either a list or a single file """

    if type(infiles) is str:
        print ""
        print "Single input file"
        tmplist = [infiles,]
    elif type(infiles) is list:
        print ""
        print "Input file list with %d members" % (len(infiles))
        tmplist = infiles
    else:
        print ""
        print "Warning.  Input frames need to be either a list of files "
        print " (python type==list) or a single input file name."
        print ""
        return

    """ Get the median fits file """
    try:
        medfits = pf.open(medfits)
    except:
        print ""
        print "ERROR: Could not open %s" % medfits
        print ""
        return

    """ Loop through the input files """
    for f in infiles:
        """ Load input file data """
        try:
            infits = pf.open(f)
        except:
            print ""
            print "ERROR: Could not open %s" % f
            print ""
            return
        indat = infits[0].data
        hdr = infits[0].header

        """ Load the weight file data, which will get updated """
        whtfile = f.replace('fits','weight.fits')
        try:
            whtfits = pf.open(whtfile,mode='update')
        except:
            print ""
            print "ERROR: Could not open %s" % whtfile
            print ""
            infits.close()
            return

        """ Set up the relevant region to be examined and the flux scale"""
        x1 = hdr['comin1'] - 1
        y1 = hdr['comin2'] - 1
        x2 = x1 + data.shape[1]
        y2 = y1 + data.shape[0]
        fscal = hdr['flxscale']

        """ 
        Make the cutout of the median-stack image and then take the 
        difference between this file and the scaled input individual file
        """
        meddat = medfits[0].data[y1:y2,x1:x2].copy()
        diff = meddat - indat * fscal

        """ Get an estimate of the RMS noise in the data """
        rms = sqrt(meddat[meddat>0])
        # Need to talk to James to interpret his idl code

        """ Close the files for this loop """
        infits.close()
        whtfits.flush()
        del meddat,diff


