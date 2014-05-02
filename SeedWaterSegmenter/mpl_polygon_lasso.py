# Original code from: http://www.mail-archive.com/matplotlib-devel@lists.sourceforge.net/msg03637.html
# Original author, Eric Bruning
# Modified by David Mashburn to allow free-drawing as well as click to form polygon

# A modified version of Eric Bruning lasso example
# (just made the callback compatible with standard lasso tool)
# and fiddled with some other stuff...
import numpy

import matplotlib
matplotlib.use('WxAgg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Widget
from matplotlib.path import Path

def points_inside_poly(points,vertices):
    '''Replacement for matplotlib.nxutils.points_inside_poly'''
    return Path(vertices).contains_points(points)

class PolyLasso(Widget):
    """ 
    A lasso widget that allows the user to define a lasso region by
    clicking to form a polygon.
    """
    
    def __init__(self, ax, callback=None, 
                 line_to_next=True, 
                 useblit=True):
        """
        Create a new lasso.
        
        *ax* is the axis on which to draw
        
        *callback* is a function that will be called with arguments (*ax*, *line*, *verts*).
        *verts* is a list of (x,y) pairs that define the polygon, with the first and 
        last entries equal.
        
        *line_to_next* controls whether or not a line is drawn to the current 
        cursor position as the polygon is created
        
        *useblit* = *True* is the only thorougly tested option.
        
        """
        self.axes = ax
        self.figure = ax.figure
        self.canvas = self.figure.canvas
        self.useblit = useblit
        self.line_to_next = line_to_next
        if useblit:
            self.background = self.canvas.copy_from_bbox(self.axes.bbox)

        
        self.verts = []
        self.line = None
        self.callback = callback
        self.cids = []
        self.cids.append(self.canvas.mpl_connect('button_release_event', self.onrelease))
        self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
        self.cids.append(self.canvas.mpl_connect('draw_event', self.ondraw))
        
    def ondraw(self, event):
        """ draw_event callback, to take care of some blit artifacts """
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        if self.line:
            self.axes.draw_artist(self.line)
        self.canvas.blit(self.axes.bbox)
        
    def do_callback(self, event):
        """ idle_event callback after polygon is finalized. """
        # Clear out callbacks. 
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.callback(self.verts)
        self.cleanup()
        
    def cleanup(self):
        """ Remove the lasso line. """
        # Clear out callbacks
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.axes.lines.remove(self.line)
        self.canvas.draw()
        
    def finalize(self):
        """ Called when user makes the final right click """
        # all done with callbacks pertaining to adding verts to the polygon
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.cids = []
        
        # Greater than three verts, since a triangle
        # has a duplicate point to close the poly.
        if 1:#len(self.verts)>3:
            # Matplotlib will not draw the closed polygon until we step completely
            # out of this routine. It is desirable to see the closed polygon on screen
            # in the event that *callback* takes a long time. So, delay callback until
            # we get an idle event, which will signal that the closed polygon has also
            # been drawn. 
            # You know, if the user makes a point or line segment, let them get no points...
            self.cids.append(self.canvas.mpl_connect('idle_event', self.do_callback))
        #else:
        #    print 'Need at least three vertices to make a polygon'
        #    self.cleanup()
                    
    def draw_update(self):
        """ Adjust the polygon line, and blit it to the screen """
        self.line.set_data(zip(*self.verts))
        if self.useblit:
            self.canvas.restore_region(self.background)
            self.axes.draw_artist(self.line)
            self.canvas.blit(self.axes.bbox)
        else:
            self.canvas.draw_idle()
    
    def onmove(self, event):
        """ Update the next vertex position """
        if self.line == None:
            # Deal with the first click
            self.line=Line2D([event.xdata], [event.ydata], 
                        linestyle='-', marker=None, markeredgecolor='white', color='white', lw=2, ms=4, animated=True)
            self.axes.add_line(self.line)
            self.verts.append((event.xdata, event.ydata))
        if event.button == None: # Move line around 
            self.verts[-1] = ((event.xdata, event.ydata))
        elif event.button == 1: # Free Hand Drawing
            self.verts[-1] = (event.xdata, event.ydata)
            self.verts.append((event.xdata, event.ydata))
        if self.line_to_next:
            self.draw_update()
    
    def onrelease(self, event):
        """ User clicked the mouse. Add a vertex or finalize depending on mouse button. """
        if self.verts is None: return
        if event.inaxes != self.axes: return
        
        if event.button == 3:
            # Right click - close the polygon
            # Set the dummy point to the first point
            self.verts[-1] = self.verts[0]
            self.draw_update()
            self.finalize()
            return

        # The rest pertains to left click only
        if event.button != 1: return
        
        if self.line == None:
            # Deal with the first click
            self.line=Line2D([event.xdata], [event.ydata], 
                        linestyle='-', marker=None, markeredgecolor='white', color='white', lw=2, ms=4, animated=True)
            self.axes.add_line(self.line)
            self.verts.append((event.xdata, event.ydata))
            
        # finalize vertex at this click, set up a new on that changes as mouse moves
        self.verts[-1] = (event.xdata, event.ydata)
        self.verts.append((event.xdata, event.ydata))

        self.draw_update()

class manager(object):
    def __init__(self):
        self.x = numpy.arange(10.0)
        self.y = self.x**2.0
        self.charge = numpy.zeros_like(self.x)
        
        self.f = plt.figure()
        ax = self.f.add_subplot(111)
        ax.set_axis_bgcolor('gray')
        self.sc = ax.scatter(self.x, self.y, c=self.charge, vmin=-1, vmax=1, edgecolor='none')
        # self.sc.set_clim(-1, 1)
        
        self.lasso = PolyLasso(ax, self.callback)
        self.f.canvas.widgetlock(self.lasso)
    
    def callback(self, verts):
        
        #print verts
        mask = points_inside_poly(zip(self.x, self.y), verts) == 1
        #print self.x[mask]
        #print self.y[mask]
        self.charge[mask] = -1
        #print self.charge
        #self.lasso_line = lasso_line
        
        # not actually necessary ... scatter stores a ref to charge array
        # self.sc.set_array(self.charge)
        
        self.f.canvas.widgetlock.release(self.lasso)
    
    
if __name__ == '__main__':
    import wx
    
    plt.ioff()
    m = manager()
    plt.show()
    
    app=wx.App();wx.MessageBox('Close to finish')
