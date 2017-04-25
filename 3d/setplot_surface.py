
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
import os, shutil
from matplotlib import rcParams
from mapping import Mapping
import dtopotools_horiz_okada_and_1d as dtopotools
reload(dtopotools)
from clawpack.geoclaw.data import LAT2METER


length_scale = 1.0e-3 # m to km
xlimits = [-150.e3*length_scale,200.e3*length_scale]
ylimits = [-87.5e3*length_scale,87.5e3*length_scale]
rcParams['contour.negative_linestyle'] = 'solid'

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """
    slice_number = 1 # set to surface slice number
    tmpdir = os.path.abspath(os.curdir)
    os.chdir(plotdata.outdir)
    for filename in os.listdir('.'):
        if (filename.startswith('slice_%d' % slice_number)):
            shutil.copyfile(filename,filename.replace('slice_%d' % slice_number,'fort',1))

    fault = dtopotools.Fault()
    fault.read('fault.data')

    os.chdir(tmpdir)

    from clawpack.visclaw import colormaps

    xc = np.linspace(-150e3,200e3,350)
    yc = np.linspace(-87.5e3,87.5e3,175)
    fault.create_dtopography(xc/LAT2METER,yc/LAT2METER,[1.0])
    Xc,Yc = np.meshgrid(xc*length_scale,yc*length_scale)

    okada_max = np.max(fault.dtopo.dZ[0,:,:])
    clevels = np.linspace(-okada_max,okada_max,10)

    plotdata.clearfigures()  # clear any old figures,axes,items data
    #plotdata.format = 'binary'

    def mapc2p(xc,yc):
        return xc*length_scale,yc*length_scale

    def plot_okada_contour(current_data):
        from pylab  import contour, ylabel, tick_params

        contour(Xc,Yc,fault.dtopo.dZ[0,:,:],levels=clevels,colors='r',linewidths=2)
        tick_params(labelsize=25)
        ylabel('kilometers',fontsize=25)

    # Figure for vertical displacement
    plotfigure = plotdata.new_plotfigure(name='surface', figno=1)
    plotfigure.kwargs = {'figsize':(11,6)}

    # Set axes for numerical solution:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = ''
    plotaxes.title_with_t = False
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_okada_contour

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 9
    plotitem.pcolor_cmap = 'BrBG'
    plotitem.pcolor_cmin = clevels[0]
    plotitem.pcolor_cmax = clevels[-1]
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 9
    plotitem.contour_colors = 'k'
    plotitem.contour_levels = clevels
    plotitem.kwargs = {'linewidths':3}
    plotitem.amr_contour_show = [0,0,0,1]
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
#    plotdata.parallel = True

    return plotdata
