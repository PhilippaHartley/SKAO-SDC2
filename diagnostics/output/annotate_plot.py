# see https://scipy-cookbook.readthedocs.io/items/Matplotlib_Interactive_Plotting.html


import math
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
from matplotlib.colors import LogNorm
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib import rc
import matplotlib.ticker as ticker




class accumulator(object):
    """ Provides for event callbacks for matplotlib drag/release events and
        axis limit changes by accumulating a series of event occurrences.
        Produces a single call to func after a user interacts with the plot.

        Sample usage:

        from pylab import figure, show

        def simple(count):
            print "update ", count
        a = accumulator(simple)

        f=figure()
        ax=f.add_subplot(111)
        plt=ax.plot(range(10))
        f.canvas.mpl_connect('draw_event', a.draw_event)
        f.canvas.mpl_connect('button_release_event', a.mouse_up_event)
        f.canvas.mpl_connect('button_press_event', a.mouse_down_event)
        ax.callbacks.connect('xlim_changed', a.axis_limit_changed)
        ax.callbacks.connect('ylim_changed', a.axis_limit_changed)
        show()

        """

    def __init__(self, func):
        self.func=func
        self.reset()
        self.counter = 0
        self.mouse_up = False

    def reset(self):
        """ Reset flags after the update function is called.
            Mouse is tracked separately.
            """
        self.limits_changed = 0
        self.got_draw = False

    def axis_limit_changed(self, ax):
        self.limits_changed += 1
        self.check_status()

    def draw_event(self, event):
        self.got_draw=True
        self.check_status()

    def mouse_up_event(self, event):
        self.mouse_up = True
        self.check_status()

    def mouse_down_event(self, event):
        self.mouse_up = False

    def both_limits_changed(self):
        """ Both x and y limits changed and the mouse is up (not dragging)
            This condition takes care of the limits being reset outside of a
            dragging context, such as the view-reset (home) button on the
            Matplotlib standard toolbar.
            """
        return (self.limits_changed >= 2) & self.mouse_up

    def interaction_complete(self):
        """ x, y, or both limits changed, and the mouse is up (not dragging).
            Also checks if matplotlib has done its final redraw of the screen,
            which comes after the call to *both* set_xlim and set_ylim
            have been triggered. The check for the draw event is the crucial
            step in not producing two calls to self.func.
        """
        return (self.limits_changed>0) & self.got_draw & self.mouse_up

    def check_status(self):
        if self.both_limits_changed() | self.interaction_complete():
            self.func(self.counter)
            self.reset()
            self.counter += 1 



class AnnoteFinder(object):
    """callback for matplotlib to display an annotation when points are
    clicked on.  The point which is closest to the click and within
    xtol and ytol is identified.
    
    Register this function like this:
    
    scatter(xdata, ydata)
    af = AnnoteFinder(xdata, ydata, annotes)
    connect('button_press_event', af)
    """

    def __init__(self, xdata, ydata, annotes, annotes_names, annotes_units, ax=None, xtol=None, ytol=None):
        self.data = list(zip(xdata, ydata, annotes))
        if xtol is None:
            xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
        if ytol is None:
            ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
        self.xtol = xtol
        self.ytol = ytol
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.drawnAnnotations = {}
        self.links = []

    def distance(self, x1, x2, y1, y2):
        """
        return the distance between two points
        """
        return(math.sqrt((x1 - x2)**2 + (y1 - y2)**2))

    def __call__(self, event):

    

        tb = plt.get_current_fig_manager().toolbar
        if event.button==1 and event.inaxes and tb.mode == '':

            if event.inaxes:

                clickX = event.xdata
                clickY = event.ydata
                if (self.ax is None) or (self.ax is event.inaxes):
                    annotes = []
                    # print(event.xdata, event.ydata)
                    for x, y, a in self.data:
                        # print(x, y, a)
                        if ((clickX-self.xtol < x < clickX+self.xtol) and
                                (clickY-self.ytol < y < clickY+self.ytol)):
                            annotes.append(
                                (self.distance(x, clickX, y, clickY), x, y, a))
                    if annotes:
                        annotes.sort()
                       
                        distance, x, y, annote = annotes[0]
                        self.drawAnnote(event.inaxes, x, y, annote, annotes_names, annotes_units)
                        for l in self.links:
                            print ('l:',l)
                            l.drawSpecificAnnote(annote)






    def drawAnnote(self, ax, x, y, annote, annotes_names, annotes_units):
        """
        Draw the annotation on the plot
        """



        xmin, xmax, ymin, ymax = ax.axis()
        plotsize_y = (ymax)-((ymax-ymin)*0.1)
        y_interval = ((ymax-ymin)*0.9)/10
        for key in self.drawnAnnotations:
            markers = self.drawnAnnotations[key]
            for m in markers:             
                m.set_visible(False)

            self.ax.figure.canvas.draw_idle()

        else:

            xpos = xmax+((xmax-xmin)*0.05)

            self.drawnAnnotations[(x, y)] = []
            for annotation in range(len(annote)):
                ypos = plotsize_y
                t = ax.text(xpos, ypos, annotes_names[annotation]+' '+annotes_units[annotation]+': '+annote[annotation])
                self.drawnAnnotations[(x, y)].append(t)
                plotsize_y-=y_interval
            m = ax.scatter([x], [y], s= 100, marker='+',c='red',zorder=100)
            self.drawnAnnotations[(x, y)].append(m)
            self.ax.figure.canvas.draw_idle()
            

    def drawSpecificAnnote(self, annote):
        annotesToDraw = [(x, y, a) for x, y, a in self.data if a == annote]
        for x, y, a in annotesToDraw:
            self.drawAnnote(self.ax, x, y, a)



if __name__=='__main__':



    if not len(sys.argv)==2:
        print ('Please supply the prefix of the image and truth catalogue')
        print ('Example usage for files my_prefix_image.fits and my_prefix_truthcat.fits:')
        print ('python annotate.py my_prefix')
        exit(0)

    prefix = sys.argv[1]

    image_name = prefix+'_image.fits'
    cat_name = prefix+'_truthcat.fits'

    sky = fits.getdata(image_name)
    sky_header = fits.getheader(image_name)
    RA_ref = sky_header['CRVAL1']
    Dec_ref = sky_header['CRVAL2']
    RA_delt = sky_header['CDELT1']*3600 # arcsec
    Dec_delt = sky_header['CDELT2']*3600 # arcsec

 

    cat_fits = fits.open(cat_name)

    cat = cat_fits[1].data

    # get the pixel size from header

    x = -1*cat['ra_offset']*60*60*(1/np.abs(RA_delt))
    y = cat['dec_offset']*60*60*(1/Dec_delt)
    x*= np.cos((30.197/360)*np.pi*2)

    x += sky.shape[0]*0.5 - 1#python indexing
    y += sky.shape[1]*0.5 - 1



    
    # properties to display
    RA = np.round(cat['RA'],6)
    Dec = np.round(cat['Dec'],6)
    z = cat['z']
    z = np.round(z, 3)
    i = (cat['i'])
    i = np.round(i, 3)
    DHI_arcsec=np.array([cat['new_HI_size_arcsec']]).astype(np.float)
    DHI_arcsec = np.round(DHI_arcsec,3)
    DHI_trecs_arcsec = np.array([cat['HI_size']]).astype(np.float)
    DHI_trecs_arcsec = np.round(DHI_trecs_arcsec,3)
    MHI = np.array([cat['MHI']]).astype(np.float)
    MHI = np.round(MHI, 3)
    fluxtot = np.array([cat['line_flux_integral']]).astype(np.float)
    fluxtot = np.round(fluxtot, 6)
    fluxtot_anna = np.array([cat['HI_flux']]).astype(np.float)
    fluxtot_anna = np.round(fluxtot, 6)
    PA = cat['PA']
 

    annotes= np.vstack([cat['Atlas_source'],RA, Dec, MHI, DHI_arcsec,DHI_trecs_arcsec ,z, i, fluxtot, fluxtot_anna,PA]).T
    annotes_names = ('Atlas source','RA','Dec', r'log$_{10}M_{\rm HI}$', r'$D_{\rm HI}$', r'$D_{\rm HI}$ trecs', r'$z$', 'inclination', r'total line flux', r'total line flux2', r'PA')
    annotes_units = ('','degrees', 'degrees', r'(M$_{\odot})$', '(arcsec)',  '(arcsec)','', '(degrees)', '(Jy-Hz)','(Jy-Hz)', '(degrees)')
    np.putmask(sky, sky==0, sky[sky>0].min())

    rms = 1e-6

    fig, ax = plt.subplots(figsize=(10,7))
    ax.imshow((sky), cmap = 'jet', norm=LogNorm(vmin=rms, vmax=0.5*np.max(sky)), origin = 'lower')
    ax.scatter(x,y, c = z, s = DHI_arcsec*5 ,cmap = 'RdBu_r', alpha = 0)
    ax.set_xlabel('RA (arcsec)')
    ax.set_ylabel('Dec (arcsec)')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x * 2.0)))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: ('%g') % (y * 2.0)))
    #cax = plt.axes([0.02, 0.1, 0.075, 0.8])

    pos = ax.get_position()
    pos.x0=0.1
    ax.set_position(pos)
    #fig.colorbar(im)
    af =  AnnoteFinder(x,y, annotes, annotes_names, annotes_units, ax=ax, xtol = 100, ytol = 100)
    def simple(count):
        pass
    a = accumulator(simple)
    fig.canvas.mpl_connect('button_press_event', af)
    fig.canvas.mpl_connect('draw_event', a.draw_event)
    fig.canvas.mpl_connect('button_release_event', a.mouse_up_event)
    fig.canvas.mpl_connect('button_press_event', a.mouse_down_event)
    ax.callbacks.connect('xlim_changed', a.axis_limit_changed)
    ax.callbacks.connect('ylim_changed', a.axis_limit_changed)
    plt.show()
