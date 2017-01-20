
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
import os, shutil
from mapping import Mapping
import dtopotools_horiz_okada_and_1d as dtopotools


length_scale = 1.0e-3 # m to km
xlimits = [-150.0e3*length_scale,200.0e3*length_scale]
zlimits = [-175e3*length_scale,0.0]

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """
    slice_number = 3
    os.chdir(plotdata.outdir)
    for filename in os.listdir('.'):
        if (filename.startswith('slice_%d' % slice_number)):
            shutil.copyfile(filename,filename.replace('slice_%d' % slice_number,'fort',1))

    os.chdir('..')

    fault = dtopotools.Fault()
    fault.read('fault.data')

    mapping = Mapping(fault)

    xp1 = mapping.xp1*length_scale
    xp2 = mapping.xp2*length_scale
    zp1 = mapping.zp1*length_scale
    zp2 = mapping.zp2*length_scale
    xcenter = mapping.xcenter
    ycenter = mapping.ycenter

    def mapc2p(xc,yc):
        xp,yp = mapping.mapc2p_xz(xc,yc)
        return xp*length_scale,yp*length_scale

    def plot_fault(current_data):
        from pylab import linspace, plot, tick_params
        xl = linspace(xp1,xp2,100)
        zl = linspace(zp1,zp2,100)
        plot(xl,zl,'g',linewidth=3)
        tick_params(labelsize=20)

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    #plotdata.format = 'binary'

    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:] + q[2,:,:])

    # Figure for trace(sigma)
    plotfigure = plotdata.new_plotfigure(name='fault', figno=1)
    plotfigure.kwargs = {'figsize':(10,4)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = zlimits
    plotaxes.title_with_t = False
    plotaxes.title = ''
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e6
    plotitem.pcolor_cmax = 1e6
    plotitem.add_colorbar = False
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
