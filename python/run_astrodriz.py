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
        infiles.sort()
        ifiles = infiles
        for f in ifiles:
            tmpfile = os.path.join(indir, f)
            os.system('cp -v %s %s/.' % (tmpfile, drizdir))
    elif isinstance(infiles, str):
        ifiles = glob.glob('%s/%s' % (indir, infiles))
        for f in ifiles:
            os.system('cp -v %s %s/.' % (f, drizdir))

    """ Return the number of input files """
    n_exp = len(ifiles)
    return n_exp

# ---------------------------------------------------------------------------

def set_defaults(inst, n_exp):
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
    if n_exp > 4:
        defpars['combine_type'] = 'median'
        if n_exp % 2 == 0:
            defpars['combine_nhigh'] = 1
        else:
            defpars['combine_nhigh'] = 0
        if n_exp > 7:
            defpars['combine_nhigh'] += 2
        if n_exp > 11:
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
    NOTE: the default final scales and pixfrac values are instrument-dependent
     and are set below
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
    if inst.lower() == 'acs':
        defpars['final_scale'] = 0.05
        defpars['final_pixfrac'] = 1.
    elif inst.lower() == 'wfc3uv' or inst.lower() == 'wfc3uvis':
        defpars['final_scale'] = 0.04
        defpars['final_pixfrac'] = 1.
    elif inst.lower() == 'wfc3ir':
        defpars['final_scale'] = 0.06
        defpars['final_pixfrac'] = 0.7
    else:
        raise ValueError('Instrument parameters not yet implement for %s'
                         % inst)

    """ Return the default settings """
    return defpars

# ---------------------------------------------------------------------------

def call_astrodriz(infiles, output, drizdir, **kwargs):
    """

    Does the actual call to astrodrizzle in the desired directory

    """

    os.chdir(drizdir)
    astd.AstroDrizzle(infiles, output=output, **kwargs)
    os.system('cp %s_drc_sci.fits ../%s.fits' % (output, output))
    os.system('cp %s_drc_ctx.fits ../%s_ctx.fits' % (output, output))
    os.system('cp %s_drc_wht.fits ../%s_wht.fits' % (output, output))
    os.chdir('..')

# ---------------------------------------------------------------------------

def run(infiles, indir, outroot, inst, drizdir='finalDrz', debug=False,
        **kwargs):
    """

    Prepares for the astrodrizzle call and then calls astrodrizzle

    """

    """ Prepare for the call """
    n_exp = driz_prep(infiles, indir, drizdir)
    print(n_exp)

    """ Set up the default parameters """
    adpars = set_defaults(inst, n_exp)
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
    call_astrodriz(infiles, outroot, drizdir, **adpars)
    


