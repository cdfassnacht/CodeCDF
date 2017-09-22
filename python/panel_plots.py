"""

A collection of classes that can be used to make multi-panel plots

"""

from astropy.io import fits as pf
from astropy.table import Table

# ===========================================================================

class PanelPlotter(object):
    """

    Generic class for plotting a multipanel figure

    """

    def __init__(self, infiles, rootnames=None, ralist=None,
                 declist=None, sizelist=None):

        """ Put the input file list into a table format """
        self.tab = Table([infiles,], names=('infile',))

        """ Add the additional columns to the table if requested """
        if rootnames is not None:
            self.tab['root'] = rootnames
        print(self.tab)

    # -----------------------------------------------------------------------

    def plot_setup(self, ncols, maxrows, panelsize=2.0):
        """
        Set parameters describint the size and arrangement of the panels
        """

        """ Get the total number of rows """
        ntot = len(self.infiles)
        nrows, mod = divmod(ntot, ncols)
        if mod > 0:
            nrows += 1

        """ See if we need to go onto more than one page """
        if nrows > maxrows:
            figsize = (ncols * panelsize, maxrows * panelsize)
            if outfile:
        if outfile:
            suff = '.%s' % (outfile.split('.')[-1])
            outname = outfile.replace(suff, '_1%s' % suff)
            totfigs = int(nrows / maxrows) + 1
            lastsize = (ncols * panelsize, (nplot % maxrows) * panelsize)
        else:
            figsize=(ncols * panelsize, nplot * panelsize)
            totfigs = 1
            if outfile:
                outname = outfile

        xsize = 1.0 / ncols - 0.005
        ny = min(nrows, maxrows)
        ysize = 1.0 / ny - 0.005

    # -----------------------------------------------------------------------

    def make_plot(self, ncols, maxrows, outfile=None):

        ycount = 1
        nfig = 1

        plt.figure(figsize=figsize)
        print len(self.tab)
        for i in range(nrows):

            """ 
            Check to see if we have reached the maximum number of rows for
            the output figure
            """
            if ycount == (maxrows+1):
                ycount = 1
                nfig += 1
                print ''
                if outfile:
                    plt.savefig(outname)
                    print 'Saving figure to %s' % outname
                    print ''
                    outname = outfile.replace(suff,'_%d%s'%(nfig,suff))
                else:
                    plt.show()
                print('')
                print('------------------------------------------------')
                print('')
                if nfig == totfigs:
                    plt.figure(figsize=lastsize)
                    ny = nplot % maxrows
                    ysize = 1.0 / ny - 0.005
                else:
                    plt.figure(figsize=figsize)

            """ Loop through the images """
            print xsize, ysize
            print ''
            for j in range(ncols):

                """ Define the image location """
                xstart = 0.0025 + 1.0 * (j % ncols) / ncols
                ystart = 1.0025 - 1.0 * ycount / ny
                if (j % ncols) == ncols-1:
                    ycount += 1

                """
                If there is no image for this panel, then skip this step
                """
                if band is None or band == 'None':
                    imax = plt.axes([xstart, ystart, xsize, ysize])
                    plt.text(0.5, 0.5, 'No Data', horizontalalignment='center',
                             color='k', transform=imax.transAxes, fontsize=10)
                    plt.axis('off')
                    continue
   
                """ Open the image """
                try:
                    f = imf.Image(infile, verbose=False)
                except:
                    print ''
                    print 'ERROR: Could not open %s' % infile
                    print ''
                    exit()
   
                """ Make the png images """
                imdef = 'radec'
                imsizen = lensinfo['imsize'][i]
                barlen = 1.0
                print ''
                imax = plt.axes([xstart,ystart,xsize,ysize])
                plt.axis('off')
                f.display(cmap='gaia',fmax=fmax,subimdef=imdef,dispunits='radec',
                          subimcent=[ra[i],dec[i]],
                          subimsize=[imsizen,imsizen], verbose=False)
   
   
                """ Put 1.0 arcsec scale bars on images """
                xmin,xmax = plt.xlim()
                ymin,ymax = plt.ylim()
                x1 = xmin + 0.1*(xmax-xmin)
                y1 = ymin + 0.1*(ymax-ymin)
                x2 = x1 - barlen # Since RA runs in opposite direction from x
                plt.plot([x1,x2],[y1,y1],color='w',lw=1)
                plt.xlim(xmin,xmax)
                plt.ylim(ymin,ymax)

                """ Label the plots """
                bandlab = '%s/%s' % (inst.upper(),band.upper())
                plt.text(0.05,0.9,root,horizontalalignment='left',
                         color='w',transform=imax.transAxes,fontsize=10)
                plt.text(0.95,0.9,bandlab,horizontalalignment='right',
                         color='w',transform=imax.transAxes,fontsize=10)
   
                plt.axis('off')
                del f

        print ''
        if outfile:
            plt.savefig(outname)
            print 'Saving figure to %s' % outname
            print ''
        else:
            plt.show()
