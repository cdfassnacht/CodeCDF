import numpy
import pyfits
import glob
import os
import shutil
import logging

path = "C:\Users\Astronomy Student\Desktop\m27 example data"
lum_files = glob.glob(os.path.join(path,"*.LUMINANCE.FIT"))
red_files = glob.glob(os.path.join(path,"*.RED.FIT"))
green_files = glob.glob(os.path.join(path,"*.GREEN.FIT"))
blue_files = glob.glob(os.path.join(path,"*.BLUE.FIT"))
dark_1x1_files = glob.glob(os.path.join(path,"*dark1x1*.FIT"))
dark_2x2_files = glob.glob(os.path.join(path,"*dark2x2*.FIT"))

def median_combine(input_files,output_file):
   data = []
   logging.info("Medium combining:")
   for filename in input_files:
      logging.info("...%s" % filename)
      data.append(pyfits.getdata(filename))
   stack = numpy.array(data)
   median = numpy.median(stack,axis=0)
   hdu = pyfits.PrimaryHDU(median)
   hdu.writeto(output_file)
   logging.info("Medium combine done.")

def subtract_images(a,b,output):
   logging.info("Subtracting images: '%s' - '%s' = '%s'" % (a,b,output))
   hdu = pyfits.PrimaryHDU(pyfits.getdata(a)+pyfits.getdata(b))
   hdu.writeto(output)

def align(old_files,output_path):
   # ack! can't do this with numpy! just copy for now...
   new_files = []
   logging.info("Aligning files (but not really):")
   for old_file in old_files:
      logging.info("...%s" % old_file)
      folder,filename = os.path.split(old_file)
      new_file = os.path.join(output_path,filename)
      shutil.copy(old_file,new_file)
      new_files.append(new_file)
   logging.info("Alignment done.")
   return new_files

def dark_files(old_files,new_path,dark_file):
   new_files = []
   for old_file in old_files:
      old_path,old_name = os.path.split(old_file)
      new_file = os.path.join(new_path,old_name)
      subtract_images(old_file,dark_file,new_file)
      new_files.append(new_file)
   return new_files

def main():
   os.mkdir(os.path.join(path,"work"))
   dark_1x1_master = os.path.join(path,"work","dark-1x1-master.fits")
   dark_2x2_master = os.path.join(path,"work","dark-2x2-master.fits")
   median_combine(dark_1x1_files,dark_1x1_master)
   median_combine(dark_2x2_files,dark_2x2_master)
   os.mkdir(os.path.join(path,"work","darked"))
   darked_lum_files = dark_files(lum_files,os.path.join(path,"work","darked"),
                                 dark_1x1_master)
   darked_red_files = dark_files(red_files,os.path.join(path,"work","darked"),
                                 dark_2x2_master)
   darked_green_files = dark_files(green_files,os.path.join(path,"work","darked"),
                                   dark_2x2_master)
   darked_blue_files = dark_files(blue_files,os.path.join(path,"work","darked"),
                                  dark_2x2_master)
   os.mkdir(os.path.join(path,"work","aligned"))
   aligned_lum_files = align(darked_lum_files,os.path.join(path,"work","aligned"))
   aligned_red_files = align(darked_red_files,os.path.join(path,"work","aligned"))
   aligned_blue_files = align(darked_blue_files,os.path.join(path,"work","aligned"))
   aligned_green_files = align(darked_green_files,os.path.join(path,"work","aligned"))
   median_combine(aligned_lum_files,os.path.join(path,"work","luminance.fits"))
   median_combine(aligned_red_files,os.path.join(path,"work","red.fits"))
   median_combine(aligned_blue_files,os.path.join(path,"work","blue.fits"))
   median_combine(aligned_green_files,os.path.join(path,"work","green.fits"))

if __name__ == "__main__":
   main()
