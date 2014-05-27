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
from math import sqrt

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

def make_wht_for_final(infiles, medfile, nsig, inwht_suff='.weight.fits',
                       outwht_suff='_wht.fits', flag_posonly=False, 
                       medwhtfile='default'):
    """
    Creates a weight file for each input image to a final swarp call.
    This weight file will assign zero weight to pixels that differ by
    more than nsig*rms from the median-stacked image that was created
    in an earlier step.
    Therefore, this function is a lot like the "blot" step in multidrizzle.

    Inputs:
       infiles      - list of individual input fits files that will be compared
                       to the median stack.  These will probably be called 
                       *resamp.fits
       medfile      - median-stacked fits file
       nsig         - minimum sigma difference between an individual input image
                      and the median stacked image that will cause a pixel to be
                      flagged.
       flag_posonly - Sets flagging behavior.  If flag_posonly=True then only
                      flag pixels where the value in the individual image is
                      more than nsig*rms _larger_ than the value in the
                      median-stacked image.
                      If False, then flag pixels for which the absolute
                      value of the difference is more than nsig*rms
                      Explanation: Set to True to flag, e.g., cosmic rays in
                      the individual image and to avoid flagging
       medwhtfile   - Name of the weight file associated with medfile.  The
                      default value ('default') will take the name of the
                      median-stacked image (medfile) and replace '.fits' with
                      '_wht.fits'
    """

    import imfuncs as imf

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

    """ Set up the weight file associated with the median-stacked image """
    if medwhtfile == 'default':
        medwhtfile = medfile.replace('.fits','_wht.fits')

    """ First loop through the input files to get information """
    gain = n.zeros(len(tmplist))
    bkgd = n.zeros(len(tmplist))
    fscal = n.zeros(len(tmplist))
    x1 = n.zeros(len(tmplist))
    y1 = n.zeros(len(tmplist))
    rmssky = n.zeros(len(tmplist))

    print ''
    print ' Input file                  gain  <bkgd>   fscal  rms_sky rms_scal'
    print '---------------------------- ---- -------- ------- ------- --------'
    for i in range(len(tmplist)):

        """ Load the header information """
        f = tmplist[i]
        try:
            hdr = pf.getheader(f)
        except:
            print ""
            print "ERROR: Could not open %s" % f
            print ""
            return

        origfile = f.replace('resamp.fits','fits')
        try:
            orighdr = pf.getheader(origfile)
        except:
            print ""
            print "ERROR: Could not open %s" % origfile
            print ""
            return

        """ Set the relevant values """
        x1[i]     = hdr['comin1'] - 1
        y1[i]     = hdr['comin2'] - 1
        bkgd[i]   = hdr['backmean']
        fscal[i]  = hdr['flxscale']
        gain[i]   = orighdr['gain']
        rmssky[i] = sqrt(bkgd[i] / gain[i])

        """ Print out the relevant information and clean up """
        print '%-28s %4.2f %8.2f %7.5f %7.3f %8.3f' \
            % (f[:-5],gain[i],bkgd[i],fscal[i],rmssky[i],rmssky[i]*fscal[i])
        del hdr,orighdr

    gainmean = gain.mean()
    print '-------------------------------------------------------------------'
    print ' %3d input files             gain  <bkgd>   fscal  rms_sky rms_scal'\
        % (len(tmplist))
    print '---------------------------- ---- -------- ------- ------- --------'
    print '    Mean values              %4.2f %8.2f %7.5f %7.3f %8.3f' % \
        (gainmean,bkgd.mean(),fscal.mean(),rmssky.mean(), \
             rmssky.mean()*fscal.mean())

    """ Open up the median file and its associated weight image"""

    """ Loop through the input files """
    epsil = 1.e-20
    print ''
    print 'Updating weight files, using diff > %d sigma' % nsig
    print '--------------------------------------------------------------------'
    for i in range(len(tmplist)):
        """ Set up file names """
        f = tmplist[i]
        inwhtfile = f.replace('.fits',inwht_suff)
        outwhtfile = f.replace('.fits',outwht_suff)

        """ Load input resamp.fits file """
        try:
            indat,hdr = pf.getdata(f,header=True)
        except:
            print ""
            print "ERROR: Could not open %s" % f
            print ""
            return

        """ Load the associated weight file data """
        try:
            inwht,whthdr = pf.getdata(inwhtfile,header=True)
        except:
            print ""
            print "ERROR: Could not open %s" % inwhtfile
            print ""
            return


        """ Set up the relevant region to be examined """
        x2 = x1[i] + indat.shape[1]
        y2 = y1[i] + indat.shape[0]

        """ 
        Make the cutout of the median-stack image and then take the 
        difference between this file and the scaled input individual file
        """
        meddat = (pf.getdata(medfile))[y1[i]:y2,x1[i]:x2]
        medwht = (pf.getdata(medwhtfile))[y1[i]:y2,x1[i]:x2]
        diff = n.zeros(indat.shape)
        whtmask = inwht>0
        diff[whtmask] = indat[whtmask] * fscal[i] - meddat[whtmask]
        del whtmask

        """ 
        Get an estimate of the RMS noise in the data.  Since we are creating
        a difference, indat - meddat, the final rms will be
            sigma_diff^2 = sigma_ind^2 + sigma_med^2
        The input variances will come from:

            sigma_ind^2 = (inddat+bkgd) / gain
              - inddat is a background-subtracted image with units of ADU.
                Therefore, add back in the background and then divide by
                the gain to get the variance in each pixel.
                This is because sigma_ADU^2 = ADU/gain, assuming that
                N_e = gain * ADU is a Poisson process.

            sigma_med^2 = (1. / medwht) + [for SNR>1: meddat/gain] 
              - Swarp produces an inverse-variance weight map
                NOTE: However, this does NOT include Poisson noise from the
                objects in the image.  This is because computes this Poisson
                noise when it does its object detections, and therefore 
                expects any input weight file to not include the Poisson noise.

        """
        medvar = n.zeros(indat.shape)
        medmask = (inwht>0) & (n.absolute(medwht>epsil))
        medvar[medmask] = 1. / medwht[medmask]
        del medmask
        medmask = meddat > n.sqrt(medvar)
        medvar[medmask] += meddat[medmask] / gainmean
        del medmask,medwht,meddat

        indvar = n.zeros(indat.shape)
        mask = (inwht>0) & ((indat + bkgd[i]) > 0.)
        indvar[mask] = fscal[i]**2 * (indat[mask]+bkgd[i])/gain[i]
        rms = n.sqrt(medvar + indvar)
        del medvar,indvar,mask,indat

        """ 
        Flag pixels that deviate by more than nsig sigma from the 
        median-stacked image
        """
        if flag_posonly:
            blotmask = diff > nsig*rms
        else:
            blotmask = n.absolute(diff) > nsig*rms
        print '%-38s %9d' \
            % (f[:-5],blotmask.sum())

        """ Debugging step(s) """
        #foodat = n.ones(diff.shape)
        #foodat[blotmask] = 0
        #fooname = f.replace('.fits','_foo.fits')
        #pf.PrimaryHDU(foodat).writeto(fooname,clobber=True)
        #del foodat

        #rmsname = f.replace('.fits','_rms.fits')
        #pf.PrimaryHDU(rms,hdr).writeto(rmsname,clobber=True)

        snr = n.zeros(diff.shape)
        rmsmask = rms>0
        snr[rmsmask] = diff[rmsmask] / rms[rmsmask]
        snrname = f.replace('.fits','_snr.fits')
        pf.PrimaryHDU(snr,hdr).writeto(snrname,clobber=True)
        del snr

        """ Write out a new weight file with the newly-flagged pixels """
        inwht[blotmask] = 0
        print '%s -> %s' % (inwhtfile,outwhtfile)
        pf.PrimaryHDU(inwht,whthdr).writeto(outwhtfile,clobber=True)

        """ Close the files for this loop """
        del diff,blotmask,rms,inwht


