import numpy as np
import math 
import cmath
import os
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from scipy.interpolate import CloughTocher2DInterpolator
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

import matplotlib.ticker as ticker
plt.switch_backend("agg")


matplotlib.rc('font', family='STIXGeneral')
from fractions import Fraction
from matplotlib.ticker import ScalarFormatter
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format


top = cm.get_cmap('Oranges_r', 128) # r means reversed version
bottom = cm.get_cmap('Blues', 128)# combine it all
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
orange_blue = ListedColormap(newcolors, name='OrangeBlue')

def read_pos(nat, nat_M_top, nat_M_bot, nat_X):
    f1=open("POSCAR")
    ave_cen0 = np.zeros(3, dtype=float)
    ave_cen1 = np.zeros(3, dtype=float)
    nat_atoms = np.zeros(2, dtype=int)

    adir0 = np.mat(np.zeros((3,3), dtype=float))
    adir1=np.mat(np.zeros((3,3), dtype=float))
    pos0 = np.mat(np.zeros((nat,3), dtype=float))
    cpos0 = np.mat(np.zeros((nat,3), dtype=float))

    f1.readline()
    line=f1.readline()
    scale=[float(line.split()[i]) for i in range(0,1)][0]
    print("scale:",scale)
    for j in range(0, 3):
        line=f1.readline()
        adir0[j,:]=[float(line.split()[i]) for i in range(0,3)]
    adir0[:,:]=scale * adir0[:,:]
    print("adir0:")
    print(adir0)
    f1.readline()
    line=f1.readline()
    nat_atoms[:] = [int(line.split()[i]) for i in range(0,2)]
    f1.readline()
    for j in range(0, nat):
        line=f1.readline()
        pos0[j,:] = [float(line.split()[i]) for i in range(0,3)]
        cpos0[j,:] = pos0[j,:] * adir0[:,:]
    f1.close()
    print("nat_atoms:", nat_atoms)
    return adir0, pos0

def read_pot(grid):
    pot_2d_top_sm = np.zeros( (grid[0], grid[1]), dtype=float)
    pot_2d_bot_sm = np.zeros((grid[0], grid[1]), dtype=float)

    f1=open("pot_2d_top_sm.dat", "r")
    f2=open("pot_2d_bot_sm.dat", "r")
    for ix in range(0, grid[0]):
        for iy in range(0, grid[1]):
            line=f1.readline()
            pot_2d_top_sm[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
            line=f2.readline()
            pot_2d_bot_sm[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
    f1.close()
    f2.close()
    return pot_2d_top_sm, pot_2d_bot_sm

def plot_pot(adir0, grid, pot_2d, sflag):
    adir0 = adir0 / adir0[0, 0]

    xtop = np.linspace(0, 1, grid[0])
    ytop = np.linspace(0, 1, grid[1])
    pot_diag = np.zeros((grid[0], 2), dtype=float)
    vtmp = np.zeros(2, dtype=float)

    for i1 in range(0, grid[0]):
        vtmp[:] = [xtop[i1], xtop[i1]] * adir0[0:2, 0:2]
        pot_diag[i1, 0] = np.sum(np.square(vtmp))**0.5 
        pot_diag[i1, 1] = pot_2d[i1, i1]


    pot_2d_sup = np.zeros((6 * grid[0], 6*grid[1]), dtype=float)
    print("shape:", pot_2d_sup.shape)
    print("shape pot_2d:", pot_2d.shape)


    for i1 in range(-3, 3):
        for i2 in range(-3, 3):
            pot_2d_sup[(i1+3) * grid[0]:(i1+4)*grid[0], (i2+3)*grid[1]:(i2+4)*grid[1] ] = pot_2d[:, :]

    xtop = np.linspace(-3, 3, 6*grid[0])
    ytop = np.linspace(-3, 3, 6*grid[1])
    X, Y = np.meshgrid(xtop, ytop)

    X = X + 0.5 * Y 
    Y = Y *3**0.5/2.0 

    X = X  * adir0[0,0]
    Y = Y * adir0[0, 0]


    fig, ax1 = plt.subplots()
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params( width=1.5, direction='in', length=8, pad=8)  

    nfont=20
#ax1.tick_params(axis='x', which='minor', direction='in', length=4, labelsize=nfont)    
    #fig.set_size_inches(9,6.75)
    #ax1.yaxis.major.formatter._useMathText = True

#for j in range(0, nbnd):
#for j in range(0, nbnd):
    nfont=20
    plt.rcParams['mathtext.fontset'] = 'stix'         

    #ax1.margins(y=.1, x=.1)

    plt.gca().set_aspect('equal', adjustable='box')
#plt.xlabel("x")
#plt.ylabel("y")
    plt.axis('off')

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 5)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#292C76", "#E6E7E8", "#EF4136"])
    cset = plt.contourf(X,Y, pot_2d_sup.T * 1000, 100, cmap=cmap, vmin=-18.7, vmax=18.7 ) # vmin=0, vmax=21) #cmap="RdBu_r",  vmin=0, vmax=21) #cmap = orange_blue)
    #cset = plt.contourf(X,Y, pot_2d_sup.T * 1000, 100, cmap="coolwarm", alpha=0.8, antialiased=True) # vmin=-19, vmax=19) #cmap = orange_blue)
    #cset = plt.contourf(X,Y, pot_2d_sup.T * 1000, 100, cmap="coolwarm", alpha=0.8, antialiased=True) # vmin=-19, vmax=19) #cmap = orange_blue)
    plt.gca().set_aspect('equal', adjustable='box')
#contour = plt.contour(Xtop,Ytop,ztop,8,colors='k')
#plt.clabel(contour, fontsize=10, colors='k')
    plt.xticks()
    plt.yticks()
    nfont=20
    #plt.xlim(-adir0[0, 0], adir0[0, 0]*1)
    #plt.ylim(-adir0[0, 0] * 3**0.5/2*1.2, adir0[0, 0] * 3**0.5/2*1.2)
    #plt.xlim([-adir0[0, 0], -adir0[0, 0] + 3 * adir0[0, 0]])
    #plt.ylim([-adir0[0, 0]* 3**0.5, -adir0[0, 0] * 3**0.5 + 1.5 * 3**0.5 * adir0[0, 0]])


    plt.xlim([-adir0[0, 0], -adir0[0, 0] + 3 * adir0[0, 0]])
    plt.ylim([-adir0[0, 0]* 3**0.5, -adir0[0, 0] * 3**0.5 + 1.5 * 3**0.5 * adir0[0, 0]])



    orig = [ 0.5 * adir0[0, 0], 3**0.5/2 * adir0[0, 0]]
    plt.plot([0.0-orig[0], adir0[0,0]-orig[0]], [0.0-orig[1], 0.0-orig[1]], ':', color='k', linewidth=1)
    plt.plot([0.0-orig[0], adir0[1, 0]-orig[0]], [0.0-orig[1], adir0[1,1]-orig[1]], ':', color='k', linewidth=1)
    plt.plot([adir0[0,0]-orig[0], adir0[0, 0]+ adir0[1,0]-orig[0]], [0.0-orig[1], adir0[1,1]-orig[1]], ':', color='k', linewidth=1)
    plt.plot([adir0[1,0]-orig[0], adir0[0, 0]+ adir0[1,0]-orig[0]], [adir0[1,1]-orig[1], adir0[1,1]-orig[1]], ':', color='k', linewidth=1)
    
    cb=plt.colorbar(cset,fraction=0.025, pad=0.083)
    tick_locator = ticker.MaxNLocator(nbins=5)  # colorbar上的刻度值个数
    cb.locator = tick_locator
    cb.ax.yaxis.major.formatter._useMathText = True
    cb.ax.yaxis.major.formatter._useOffset = True
    #cb.ax.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f")) #OOMFormatter(-3, "%1.1f"))
    #cb.set_ticks([0.0, 1.0, 2.0, 3.0 ])
#cb.set_ticklabels([-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5, 2.0, 2.5])
    cb.update_ticks()
    #plt.clim(-19, 19)
    plt.title('$\Delta v_{H}$(r)', size=nfont+4)
    cb.ax.set_title('(meV)', size=nfont+2)
    cb.ax.tick_params(labelsize=nfont)
    plt.rcParams['font.sans-serif']='stix'

#plt.rcParams['mathtext.fontset'] = 'stix'         
#cb.set_label('Strain-top')
#plt.xlabel("x")
#plt.ylabel("y")
    #cb.remove()
    #fig.tight_layout()
    plt.axis('off')
    #plt.title(sflag)
    plt.savefig('fig-paper-pot-'+sflag+'.png',dpi=500) #  transparent=True)
    return pot_diag

def plot_pot_diag(adir0, pot_diag, sflag):
    adir0 = adir0 / adir0[0,0]
    fig, ax1 = plt.subplots()
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params( width=1.5, direction='in', length=8, pad=8)  

    nfont = 30
    ax1.tick_params(axis='x', which='minor', direction='in', length=4, labelsize=nfont)    
    fig.set_size_inches(9, 6)
    ax1.xaxis.major.formatter._useMathText = True
    plt.rcParams['mathtext.fontset'] = 'stix'         

#for j in range(0, nbnd):
    
    ax1.plot(pot_diag[:, 0], pot_diag[:, 1]*1000, color="tomato",linestyle="solid", linewidth=2) #label="on-site-nsc")
    ax1.set_ylabel("v$_{H}$(r) (meV)", fontsize=nfont+2, fontname="STIXGeneral")
    ax1.legend(frameon=False, fontsize=nfont,labelspacing = 0.0, loc='lower left', handlelength=1, handletextpad=0.3)


    plt.xticks([])
    #plt.xlim(0, adir0[0, 0] *1.5)
    plt.text(0.0, -2.5, 'MM',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")
    #plt.text(kx[20], -108, 'm',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")
    #plt.text(kx[40] , -108, '$\kappa$',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")
    #plt.text(kx[60], -108, '$\gamma$',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")
    #plt.axvline(x=kx[20], color='k', alpha=0.3, linestyle='-', linewidth=1.0)
    #plt.axvline(x=kx[40], color='k', alpha=0.3, linestyle='-', linewidth=1.0)
#plt.axhline(y=0.0, color='k', alpha=0.5, linestyle='--', linewidth=1.0)
    #plt.ylim(-100,10)
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(1.5)
    ax1.tick_params( width=1.5)
    plt.yticks(fontsize=nfont, fontname="STIXGeneral")
    plt.xticks(fontsize=nfont, fontname="STIXGeneral" ) # fontname='Arial')
    fig.tight_layout()
    plt.savefig("plot-pot-diag"+sflag+".png", dpi=500)
    return 




nat =        1324 +  662
nat_M_top = 331
nat_M_bot = 331
nat_X = 1324 

adir0, pos0 = read_pos(nat, nat_M_top, nat_M_bot, nat_X) 


grid = np.zeros(3, dtype=int)

f0=open("plot-pot.input", "r")
line=f0.readline()
grid[:]=[float(line.split()[i]) for i in range(0,3)]
f0.close()

print("grid:", grid[0], grid[1])
pot_2d_top_sm, pot_2d_bot_sm = read_pot(grid)

min_val = min(np.reshape(pot_2d_top_sm, (grid[0]*grid[1], -1)))
plot_pot(adir0, grid, pot_2d_top_sm-min_val, 'top')

min_val = min(np.reshape(pot_2d_bot_sm, (grid[0]*grid[1], -1)))
plot_pot(adir0, grid, pot_2d_bot_sm-min_val, 'bot')


pot_diag = plot_pot(adir0, grid, pot_2d_top_sm-pot_2d_bot_sm, 'top-bot')
#print("max:", max(np.reshape(pot_2d_top_sm-min_val, (grid[0]*grid[1], -1))))
print("max:", max(np.reshape(pot_2d_top_sm-pot_2d_bot_sm, (grid[0]*grid[1], -1))))
print("min:", min(np.reshape(pot_2d_top_sm-pot_2d_bot_sm, (grid[0]*grid[1], -1))))
plot_pot_diag(adir0, pot_diag, 'top-bot')


f=open("pot_diag-top.dat", 'w')
for i1 in range(0, grid[0]):
    np.savetxt(f, np.mat(pot_diag[i1, :]), fmt="%20.14E %20.14E")
f.close()
#min_val = min(np.reshape(pot_2d_bot_sm, (grid[0]*grid[1], -1)))
#plot_pot(adir0, grid, pot_2d_bot_sm-min_val, 'bot')

#print("max:", max(np.reshape(pot_2d_top_sm-pot_2d_bot_sm, (grid[0]*grid[1], -1))))
#print("min:", min(np.reshape(pot_2d_top_sm-pot_2d_bot_sm, (grid[0]*grid[1], -1))))

#plot_pot(adir0, grid, pot_2d_top_sm-pot_2d_bot_sm, 'top-bot')
