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

