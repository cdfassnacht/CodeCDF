import os
import sys
import glob
from drizzlepac import tweakreg
from drizzlepac import astrodrizzle as astd

# indir = sys.argv[1]
# outroot = sys.argv[2]
# 
# doCTE = False
# if len(sys.argv) > 3:
#     if sys.argv[3] == 'doCTE':
#         doCTE = True

# ---------------------------------------------------------------------------

def driz_prep(infiles, indir, drizdir='finalDrz'):
    """

    Prepares the directory for drizzling

    """
    os.system('rm -rf %s' % drizdir)
    os.system('mkdir %s' % drizdir)

    """ Get the filenames """
    if isinstance(infiles, list):
        ifiles = infiles.sort()
        for f in ifiles:
            tmpfile = os.path.join(indir, f)
            os.system('cp -v %s %s/.' % (tmpfile, drizdir))
    elif isinstance(infiles, str):
        ifiles = glob.glob('%s/%s' % (indir, infiles))
        for f in ifiles:
            os.system('cp -v %s %s/.' % (f, drizdir))

# ---------------------------------------------------------------------------

def call_astrodriz(infiles, output, final_scale, build=False, updatewcs=False,
                   static=True, skysub=True, driz_separate=True,
                   driz_sep_kernel='turbo', driz_sep_scale=None,
                   driz_sep_pixfrac=1., driz_sep_bits=0, median=True,
                   blot=True, driz_cr=True,
                   driz_cr_corr=False, driz_combine=True, final_wcs=True,
                   final_bits=0, final_kernel='square',
                   # final_kernel='lanczos3',
                   final_pixfrac=1.,final_wt_scl='exptime',
                   final_wht_type='EXP', clean=True, final_rot=0, 
                   **kwargs):
    """

    Does the actual call to astrodrizzle

    """

    astd.AstroDrizzle(infiles, output=output, final_scale=final_scale,
                      build=build, updatewcs=updatewcs, static=static,
                      skysub=skysub, driz_separate=driz_separate,
                      driz_sep_kernel=driz_sep_kernel,
                      driz_sep_scale=driz_sep_scale,
                      driz_sep_pixfrac=driz_sep_pixfrac,
                      driz_sep_bits=driz_sep_bits, median=median,
                      blot=blot, driz_cr=driz_cr,
                      driz_cr_corr=driz_cr_corr, driz_combine=driz_combine,
                      final_wcs=final_wcs, final_bits=final_bits,
                      final_kernel=final_kernel, final_pixfrac=final_pixfrac,
                      final_wt_scl=final_wt_scl, 
                      final_wht_type=final_wht_type, clean=clean,
                      final_rot=final_rot, **kwargs)

# ---------------------------------------------------------------------------

def run(infiles, indir, outroot, inst, outscale='default',
        drizdir='finalDrz', **kwargs):
    """

    Prepares for the astrodrizzle call and then calls astrodrizzle

    """

    """ Prepare for the call """
    driz_prep(infiles, indir, drizdir)

    """ Set up instrument-dependent default values """
    if inst.lower() == 'acs':
        final_scale = 0.05
    elif inst.lower() == 'wfc3uv' or inst.lower() == 'wfc3uvis':
        final_scale = 0.04
    elif inst.lower() == 'wfc3ir':
        final_scale = 0.1283
    else:
        raise ValueError('Instrument parameters not yet implement for %s'
                         % inst)

    """ Override the defaults if requested """
    if outscale.lower() != 'default':
        if isinstance(outscale, float):
            final_scale = outscale
        else:
            raise TypeError('outscale must be a float')

    """ Call the astrodrizzle routine """
    os.chdir(drizdir)
    call_astrodriz(infiles, outroot, final_scale, **kwargs)

    os.system('cp %s_drc_sci.fits ../%s.fits'%(outroot,outroot))
    os.system('cp %s_drc_ctx.fits ../%s_ctx.fits'%(outroot,outroot))
    os.system('cp %s_drc_wht.fits ../%s_wht.fits'%(outroot,outroot))
    os.chdir('..')


