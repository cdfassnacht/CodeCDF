"""

A collection of classes that can be used to make multi-panel plots

"""

import sys
from astropy.io import fits as pf
from astropy.table import Table
from matplotlib import pyplot as plt
import imfuncs as imf

# ===========================================================================

class PanelPlotter(object):
    """

    Generic class for plotting a multipanel figure

    """

    def __init__(self, infiles, ra, dec, size, fmax, lab=None, lab2=None):

        """ Put the input file list into a table format """
        self.tab = Table([infiles, ra, dec, size, fmax],
                         names=('infile', 'ra', 'dec', 'size', 'fmax'))
        if lab is not None:
            self.tab['lab'] = lab
        if lab2 is not None:
            self.tab['lab2'] = lab2

        print(self.tab)

    # -----------------------------------------------------------------------

    def plot_setup(self, ncols, maxrows, panelsize=2.0):
        """
        Set parameters describint the size and arrangement of the panels
        """

        """ Get the total number of rows """
        ntot = len(self.tab)
        self.nrows, mod = divmod(ntot, ncols)
        if mod > 0:
            self.nrows += 1

        """ See if we need to go onto more than one page """
        if self.nrows > maxrows:
            self.figsize = (ncols * panelsize, maxrows * panelsize)
            self.totfigs = int(self.nrows / maxrows) + 1
            self.lastsize = (ncols * panelsize, 
                             (self.nrows % maxrows) * panelsize)
        else:
            self.figsize=(ncols * panelsize, self.nrows * panelsize)
            self.totfigs = 1


    # -----------------------------------------------------------------------

    def make_plot(self, ncols, maxrows, outfile=None):

        """ Do some plot setup work """
        self.plot_setup(ncols, maxrows)
        xsize = 1.0 / ncols - 0.005
        ny = min(self.nrows, maxrows)
        ysize = 1.0 / ny - 0.005

        """ Set the name of the output file """
        if outfile:
            if self.nrows > maxrows:
                suff = '.%s' % (outfile.split('.')[-1])
                outname = outfile.replace(suff, '_1%s' % suff)
            else:
                outname = outfile

        """' Set up column, row, and figure counters """
        xcount = 0
        ycount = 1
        nfig = 1

        """ Open the figure and start the loop """
        plt.figure(figsize=self.figsize)

        for i in self.tab:
            """ 
            Check to see if we have reached the maximum number of rows for
            the output figure
            """
            if ycount == (maxrows + 1):
                ycount = 1
                nfig += 1
                print ''
                if outfile:
                    plt.savefig(outname)
                    print 'Saving figure to %s' % outname
                    print ''
                    outname = outfile.replace(suff,'_%d%s' % (nfig, suff))
                else:
                    plt.show()
                print('')
                print('------------------------------------------------')
                print('')
                if nfig == self.totfigs:
                    plt.figure(figsize=self.lastsize)
                    ny = self.nrows % maxrows
                    ysize = 1.0 / ny - 0.005
                else:
                    plt.figure(figsize=self.figsize)

            """ Define the image location """
            xstart = 0.0025 + 1.0 * (xcount % ncols) / ncols
            ystart = 1.0025 - 1.0 * ycount / ny
            print xcount, ycount
            print xstart, ystart
            if (xcount % ncols) == ncols - 1:
                xcount = -1
                ycount += 1

            """
            If there is no image for this panel, then skip this step
            """
            if i['infile'] is None:
                imax = plt.axes([xstart, ystart, xsize, ysize])
                plt.text(0.5, 0.5, 'No Data', horizontalalignment='center',
                         color='k', transform=imax.transAxes, fontsize=10)
                plt.axis('off')
                continue
   
            """ Open the image """
            try:
                f = imf.Image(i['infile'], verbose=True)
            except:
                print ''
                print 'ERROR: Could not open %s' % i['infile']
                print ''
                return
   
            """ Make the png images """
            imdef = 'radec'
            imsizen = i['size']
            barlen = 1.0
            print ''
            imax = plt.axes([xstart, ystart, xsize, ysize])
            plt.axis('off')
            f.display(cmap='gaia', fmax=i['fmax'], subimdef=imdef,
                      dispunits='radec', subimcent=[i['ra'],i['dec']],
                      subimsize=[i['size'],i['size']], verbose=True)
   
   
            """ Put 1.0 arcsec scale bars on images """
            xmin, xmax = plt.xlim()
            ymin, ymax = plt.ylim()
            x1 = xmin + 0.1 * (xmax - xmin)
            y1 = ymin + 0.1 * (ymax - ymin)
            x2 = x1 - barlen # Since RA runs in opposite direction from x
            plt.plot([x1, x2], [y1, y1], color='w', lw=1)
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            """ Label the plots if requested """
            lab1 = False
            lab2 = False
            if 'lab' in self.tab.colnames:
                lab1 = True
            if 'lab2' in self.tab.colnames:
                lab2 = True
            if lab1 and lab2:
                plt.text(0.05, 0.9, i['lab'], horizontalalignment='left',
                         color='w', transform=imax.transAxes, fontsize=10)
                plt.text(0.95, 0.9, i['lab2'], horizontalalignment='right',
                         color='w', transform=imax.transAxes, fontsize=10)
            else:
                lab = None
                if lab1:
                    lab = i['lab']
                elif lab2:
                    lab = i['lab2']
                if lab:
                    plt.text(0.5, 0.9, lab, horizontalalignment='center',
                             color='w', transform=imax.transAxes, fontsize=10)
   
            plt.axis('off')
            del f
            xcount += 1

        print('')
        if outfile:
            plt.savefig(outname)
            print('Saving figure to %s' % outname)
            print('')
        else:
            plt.show()
