import scipy
import sys,os
try:
   from astropy.io import fits as pyfits
except:
   import pyfits

if len(sys.argv)>1:
   files = sys.argv[1:]
else:
   dir = os.curdir
   files = os.listdir(dir)

print \
    "#   File          nx x ny         Object       t_exp    Instrument Filter  "
print \
    "#--------------  ---------  ----------------- --------- ---------- --------"

for f in files:

   # Open the fits file
   try:
      hdulist = pyfits.open(f)
   except:
      try:
         hdulist = pyfits.open(f,ignore_missing_end=True)
      except:
         print "Unable to open file %s" % f
         continue

   if f[-5:] == ".fits":
      fname = f[:-5]
   elif f[-4:] == ".FIT":
      fname = f[:-4]
   else:
      fname = f

   # Get information that is likely to be in the PHDU
   # Perhaps most important is the instrument
   hdr = hdulist[0].header

   # Instrument
   try:
      inst = hdr['instrume'].strip()
   except:
      try:
         inst = hdr['currinst'].strip()
      except:
         inst = None

   # Set some default values for instrument-specific parameters
   usecoadd  = False
   inst2     = None
   objname   = 'object'
   texpname  = 'exptime'
   texpname2 = 'exptime'
   filtname  = 'filter'
   coaddname = 'coadds'

   # Set instrument-specific parameter value
   if inst is not None:

      """ HST Instruments """
      if inst == 'ACS':
         objname = 'targname'
         filtname = 'filter2'
      if inst == 'NICMOS':
         objname = 'targname'
         inst2   = 'aperture'
      elif inst == 'WFC3':
         objname = 'targname'
      elif inst == 'WFPC2':
         objname = 'targname'
         filtname = 'filtnam1'

      ### Keck Instruments ###
      elif inst[0:3] == 'ESI':
         inst = 'ESI'
         texpname = 'ttime'
         filtname = 'dwfilnam'
      elif inst[0:4] == 'LRIS':
         inst = 'LRIS'
         texpname = 'ttime'
         filtname = 'redfilt'
      elif inst == 'NIRC2':
         usecoadd = True
         texpname = 'itime'
      elif inst == 'NIRSPEC':
         usecoadd = True
         texpname = 'itime'
         filtname = 'filname'

      ### Subaru Instruments ###
      elif inst == 'MOIRCS':
         filtname = 'filter01'
      elif inst == 'SuprimeCam':
         filtname = 'filter01'

      ### Other Instruments ###
      elif inst == 'NIRI':
         usecoadd  = True
         texpname  = 'coaddexp'
         filtname  = 'filter1'

   # Get object name
   try:
      obj = hdr[objname].strip()
   except:
      obj = '(No_object)'

   # Get data size
   if len(hdulist) > 1:
      hext = 1
   else:
      hext = 0
   try:
      n1 = hdulist[hext].header['naxis1']
   except:
      n1 = 0
   try:
      n2 = hdulist[hext].header['naxis2']
   except:
      n2 = 0

   # Get exposure time, possibly including coadds
   try:
      exptime = hdr[texpname]
   except:
      try:
         exptime = hdr[texpname2]
      except:
         exptime = float('NaN')
   if usecoadd:
      try:
         coadd = hdr[coaddname]
      except:
         coadd = 1
      texp = '%2dx%-6.2f' % (coadd,exptime)
   else:
      texp = '%7.2f  ' % exptime

   # Get filter
   try:
      filt = hdr[filtname]
   except:
      filt = 'N/A'

   # Get instrument-specific cards

   if inst == 'NICMOS':
      try:
         aper = hdr['aperture'][0:4]
      except:
         aper = ''
      inst = '%s-%s' % (inst,aper)
   if(inst == 'NIRC2'):
      try:
         cam = hdr['camname']
      except:
         cam = 'N/A'
      cam = cam[0:4]
      inst = "%s-%s" % (inst,cam)

   # Print out final info
   print "%-15s  %4dx%-4d  %-17s %-9s %-10s %s" % \
          (fname,n1,n2,obj,texp,inst,filt)

