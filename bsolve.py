# -*- coding: utf-8 -*-

# axisymmetric (RZ) particle in cell code example
#
# see https://www.particleincell.com/2015/rz-pic/ for more info

import numpy
import pylab as pl
import math

    
def plot(ax,data,scatter=False):
    pl.sca(ax)
    pl.cla()
    cf = pl.contourf(pos_z, pos_r, numpy.transpose(data),8,alpha=.75,linewidth=1,cmap='jet')
    #cf = pl.pcolormesh(pos_z, pos_r, numpy.transpose(data))
    ax.set_yticks(pos_r)
    ax.set_xticks(pos_z)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    pl.xlim(min(pos_z),max(pos_z))
    pl.ylim(min(pos_r),max(pos_r))
    ax.grid(b=True,which='both',color='k',linestyle='-')
    ax.set_aspect('equal', adjustable='box')
 #   pl.colorbar(cf,ax=pl.gca(),orientation='horizontal',shrink=0.75, pad=0.01)

    
#---------- INITIALIZATION ----------------------------------------

pl.close('all')
pos_z = [];
pos_r = [];       
b = [];         
         
#----------- END OF MAIN LOOP ------------------------
plot(pl,b)               
#Q = pl.quiver(pos_z, pos_r, numpy.transpose(efz), numpy.transpose(efr),units='xy')
pl.draw()
      

#this will block execution until figure is closed
pl.show()

