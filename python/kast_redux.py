"""
A quick-look reduction script for Kast data
"""

import numpy as n
import imfuncs as imf
import spec_simple as ss
from matplotlib import pyplot as plt

def redux_all(infile, cam, redrange=[50,211], bluerange=[920,1216],
              muorder=3, sigorder=3, fitrange=None, apmin=-4., apmax=4.,
              wavecol=None):

    """ Read in the data """
    spec2d = imf.Image(infile)

    """ Set the trim range for the input data file"""
    if cam=='red':
        xmin = redrange[0]
        xmax = redrange[1]
    elif cam=='blue':
        xmin = bluerange[0]
        xmax = bluerange[1]
    else:
        print ''
        print "The cam variable must be either 'blue' or 'red'"
        print ''
        exit()

    """ Trim and transpose the data """
    dat2 = (spec2d.hdu[0].data[:,xmin:xmax]).transpose()
    spec2d.hdu[0].data = dat2
    hdr = spec2d.hdu[0].header
    hdr['naxis1'] = dat2.shape[1]
    hdr['naxis2'] = dat2.shape[0]
    print ''

    """ Plot the 2D spectrum """
    plt.figure(1)
    plt.subplot(211)
    spec2d.display()

    """ Plot the 2D spectrum with rough sky subtraction """
    tmpsky1d = n.median(dat2,axis=0)
    tmpsky = n.tile(tmpsky1d,(dat2.shape[0],1))
    ss2 = dat2 - tmpsky
    spec2d.hdu[3].data = ss2
    hdr3 = spec2d.hdu[3].header
    hdr3['naxis1'] = ss2.shape[1]
    hdr3['naxis2'] = ss2.shape[0]    
    plt.subplot(212)
    spec2d.found_rms = False
    spec2d.display(hext=3,sighigh=3.)
    spec2d.found_rms = False

    """ Find the trace and fit a polynomial to its location and width """
    plt.figure(2)
    pos,width = ss.find_and_trace(dat2,do_subplot=True,muorder=muorder,
                                  sigorder=sigorder,fitrange=fitrange,
                                  apmin=apmin,apmax=apmax)

    """ Extract the spectrum"""
    f,v = ss.extract_spectrum(dat2,pos,width,do_subplot=True)

    """ Plot the extracted spectrum along with the wavelength """
    if wavecol is None:
        try:
            w = spec2d.hdu[5].data['wave_opt'][1,:]
        except:
            w = spec2d.hdu[5].data['wave_opt'][0,:]
    else:
        w = spec2d.hdu[5].data['wave_opt'][wavecol,:]
    plt.figure(3)
    ss.plot_spectrum_array(w,f,v)

    """ Save the spectrum """
    outspec = infile.replace('.fits','.spec')
    ss.save_spectrum(outspec,w,f,v)
