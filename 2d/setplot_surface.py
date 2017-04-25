
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from clawpack.clawutil.data import ClawData
import dtopotools_horiz_okada_and_1d as dtopotools
reload(dtopotools)
from clawpack.geoclaw.data import LAT2METER

length_scale = 1.0e-3 # m to km
xlimits = [-150e3*length_scale,200e3*length_scale]
ylimits = [-0.3,0.6]

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    fault = dtopotools.Fault()
    fault.read(plotdata.outdir + '/fault.data')
    gdata = ClawData()
    gdata.read(plotdata.outdir + '/gauges.data',force=True)
    ngauges = gdata.ngauges
    xc = np.zeros(ngauges)
    for j in range(ngauges):
        g = plotdata.getgauge(j)
        xc[j] = g.location[0]

    fault.create_dtopography(xc/LAT2METER,np.array([0.]),[1.0],y_disp=True)

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    def plot_vertical_displacement(current_data):
        from pylab import plot,zeros,xlabel,ylabel,tick_params
        t = current_data.t

        ys = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                v = 0.5*(g.q[4,k]+g.q[4,k-1])
                ys[gaugeno] += dt*v
        plot(xc[:ngauges]*length_scale,ys,linewidth=3)
        plot(xc*length_scale,fault.dtopo.dZ[0,0,:],linestyle='--',color='r',linewidth=3)
        tick_params(labelsize=25)
        xlabel('kilometers',fontsize=25)
        ylabel('meters',fontsize=25)

    def plot_horizontal_displacement(current_data):
        from pylab import plot,zeros,xlabel,ylabel,tick_params
        t = current_data.t

        xs = zeros(ngauges)
        for gaugeno in range(ngauges):
            g = plotdata.getgauge(gaugeno)
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                u = 0.5*(g.q[3,k]+g.q[3,k-1])
                xs[gaugeno] += dt*u

        plot(xc[:ngauges]*length_scale,xs,linewidth=3)
        plot(xc*length_scale,fault.dtopo.dY[0,0,:],linestyle='--',color='r',linewidth=3)
        tick_params(labelsize=25)
        xlabel('kilometers',fontsize=25)
        ylabel('meters',fontsize=25)

    # Figure for vertical displacement
    plotfigure = plotdata.new_plotfigure(name='vertical', figno=1)
    plotfigure.kwargs = {'figsize':(11,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title_with_t = False
    plotaxes.title = ''
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_vertical_displacement

    # Figure for horizontal displacement
    plotfigure = plotdata.new_plotfigure(name='horizontal', figno=2)
    plotfigure.kwargs = {'figsize':(11,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title_with_t = False
    plotaxes.title = ''
    plotaxes.scaled = False
    plotaxes.afteraxes = plot_horizontal_displacement

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
    plotdata.parallel = True

    return plotdata
