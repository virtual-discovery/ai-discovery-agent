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

def read_pot(grid, grid2):
    Delta_t = np.zeros( (grid[0], grid[1]), dtype=float)
    #pot_2d_bot_sm = np.zeros((grid[0], grid[1]), dtype=float)

    f1=open("Delta-t.dat", "r")
    #f2=open("pot_2d_bot_sm.dat", "r")
    for ix in range(0, grid[0]):
        for iy in range(0, grid[1]):
            line=f1.readline()
            Delta_t[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
            #line=f2.readline()
            #pot_2d_bot_sm[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
    f1.close()
    #f2.close()

    Delta_b = np.zeros( (grid[0], grid[1]), dtype=float)
    #pot_2d_bot_sm = np.zeros((grid[0], grid[1]), dtype=float)

    f1=open("Delta-b.dat", "r")
    #f2=open("pot_2d_bot_sm.dat", "r")
    for ix in range(0, grid[0]):
        for iy in range(0, grid[1]):
            line=f1.readline()
            Delta_b[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
            #line=f2.readline()
            #pot_2d_bot_sm[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
    f1.close()
    #f2.close()

    Delta_T = np.zeros( (grid2[0], grid2[1]), dtype=complex)
    #pot_2d_bot_sm = np.zeros((grid[0], grid[1]), dtype=float)

    f1=open("Delta-T.dat", "r")
    #f2=open("pot_2d_bot_sm.dat", "r")
    for ix in range(0, grid2[0]):
        for iy in range(0, grid2[1]):
            line=f1.readline()
            tmp1 = [float(line.split()[i]) for i in range(0,2)][0]
            tmp2 = [float(line.split()[i]) for i in range(0,2)][1]
            Delta_T[ix, iy]=  complex(tmp1, -tmp2)
            #line=f2.readline()
            #pot_2d_bot_sm[ix, iy]=[float(line.split()[i]) for i in range(0,1)][0]
    f1.close()
    #f2.close()




    return Delta_t, Delta_b, Delta_T

def plot_pot(grid, grid2, Delta_t, Delta_b, Delta_T):

    adir0 = np.zeros((2,2), dtype=float)
    #adir0[0, 0] = 1.0  
    #adir0[1, 0] =  0.5 
    #adir0[1, 1] = 3**0.5/2 
    #adir0[0, 0] = 3**0.5/2 
    #adir0[0, 1] = 0.5 
    #adir0[1, 1] = 1.0

    #adir0[0, 0] = 3**0.5/2 
    #adir0[0, 1] = 0.5 

    #adir0[1, 0] = 3**0.5/2 
    #adir0[1, 1] = -0.5 

    adir0[0, 0] = 0.5 
    adir0[0, 1] = 3**0.5/2 

    adir0[1, 0] = -0.5 
    adir0[1, 1] = 3**0.5/2 

    #adir0 = adir0 / adir0[1, 1]

    xtop = np.linspace(0, 1, grid[0])
    ytop = np.linspace(0, 1, grid[1])
   

    xtop2 = np.linspace(0, 1, grid2[0]) 
    ytop2 = np.linspace(0, 1, grid2[1])


    Delta_z_sup = np.zeros((6 * grid[0], 6*grid[1]), dtype=float)
    Delta_T_sup_real = np.zeros((6*grid2[0], 6*grid2[1]), dtype=float) 
    Delta_T_sup_imag = np.zeros((6*grid2[0], 6*grid2[1]), dtype=float)

    xy_sup2 = np.zeros((6*grid2[0], 6*grid2[1], 2), dtype=float)

    Rot = np.zeros((2, 2), dtype=float)
    deg = math.pi / 6 
    Rot[0, :] = [math.cos(deg), -math.sin(deg)]
    Rot[1, :] = [math.sin(deg), math.cos(deg)]

    

    xy_unit2 = np.zeros((grid2[0], grid2[1], 2), dtype=float)
    vtmp = np.zeros(2, dtype=float)
    for i1 in range(0, grid2[0]):
        for i2 in range(0, grid2[1]):
            #vtmp[:] = np.mat([xtop2[i1], ytop2[i2]]) * np.mat(adir0)
            #print("vtmp:", vtmp)
            xy_unit2[i1, i2, :] = [xtop2[i1], ytop2[i2]]

    
    for i1 in range(-3, 3):
        for i2 in range(-3, 3):
            Delta_z_sup[(i1+3) * grid[0]:(i1+4)*grid[0], (i2+3)*grid[1]:(i2+4)*grid[1] ] = (Delta_b[:, :]-Delta_t)/2 
            Delta_T_sup_real[(i1+3) * grid2[0]:(i1+4)*grid2[0], (i2+3)*grid2[1]:(i2+4)*grid2[1] ] = Delta_T[:, :].real 
            Delta_T_sup_imag[(i1+3) * grid2[0]:(i1+4)*grid2[0], (i2+3)*grid2[1]:(i2+4)*grid2[1] ] = Delta_T[:, :].imag


            xy_sup2[(i1+3) * grid2[0]:(i1+4)*grid2[0], (i2+3)*grid2[1]:(i2+4)*grid2[1], 0] = xy_unit2[:, :, 0] + i1 
            xy_sup2[(i1+3) * grid2[0]:(i1+4)*grid2[0], (i2+3)*grid2[1]:(i2+4)*grid2[1], 1] = xy_unit2[:, :, 1] + i2 




    data_top = np.zeros((6*grid2[0]*6*grid2[1], 4), dtype=float)
    data_top[:, 0] = np.reshape(xy_sup2[:, :, 0], ((6*grid2[0]*6*grid2[1], -1)))[:, 0]
    data_top[:, 1] = np.reshape(xy_sup2[:, :, 1], ((6*grid2[0]*6*grid2[1], -1)))[:, 0]
    data_top[:, 2] = np.reshape(Delta_T_sup_real, ((6*grid2[0]*6*grid2[1], -1)))[:, 0]/1000/2
    data_top[:, 3] = np.reshape(Delta_T_sup_imag, ((6*grid2[0]*6*grid2[1], -1)))[:, 0]/1000/2



    vtmp = np.zeros(2, dtype=float)
    ztmp = np.zeros(2, dtype=float)
    for i1 in range(0, 6*grid2[0]*6*grid2[1]):
        vtmp[0] = data_top[i1, 0] 
        vtmp[1] = data_top[i1, 1]
        ztmp = np.mat(vtmp) * np.mat(adir0)
        #print("ztmp:", ztmp[0, 0], ztmp[0, 1])
        data_top[i1, 0] = ztmp[0, 0] 
        data_top[i1, 1]  = ztmp[0, 1] 

        vtmp[0] = data_top[i1, 2] 
        vtmp[1] = data_top[i1, 3]
        ztmp = (np.mat(Rot) * np.mat(vtmp).T).T
        print("ztmp:", ztmp[0, 0], ztmp[0, 1])
        data_top[i1, 2] = ztmp[0, 0] 
        data_top[i1, 3] = ztmp[0, 1]

    xtop = np.linspace(-3, 3, 6*grid[0])
    ytop = np.linspace(-3, 3, 6*grid[1])
    X, Y = np.meshgrid(xtop, ytop)

    #Y = 0.5 * X  + Y 
    #X = X * 3**0.5/ 2
    #X = X + 0.5 * Y 
    #Y = Y *3**0.5/2.0 

    X0, Y0 = np.meshgrid(xtop, ytop) 
    
    #X = 3**0.5/2 * X0 + 3**0.5/2* Y0 
    #Y = 0.5 * X0 -0.5 * Y0  

    X = 0.5 * X0 -0.5 * Y0 
    Y = 3**0.5/2 * X0 +3**0.5/2 * Y0 

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
    fig.set_size_inches(7.5, 5)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#637493", "#CDD4DF", "#ED9F30"])
    cset = plt.contourf(X,Y, Delta_z_sup.T, 100, cmap="magma") # vmin=0, vmax=21) #cmap="RdBu_r",  vmin=0, vmax=21) #cmap = orange_blue)
    plt.quiver(data_top[:, 0], data_top[:, 1], data_top[:, 2], data_top[:, 3], units='x', pivot='tail', scale=0.2, color="w", width=0.008, headwidth=4 ) #width=1.2, headwidth=3, color="w")
    #cset = plt.contourf(X,Y, pot_2d_sup.T * 1000, 100, cmap="coolwarm", alpha=0.8, antialiased=True) # vmin=-19, vmax=19) #cmap = orange_blue)
    #cset = plt.contourf(X,Y, pot_2d_sup.T * 1000, 100, cmap="coolwarm", alpha=0.8, antialiased=True) # vmin=-19, vmax=19) #cmap = orange_blue)
    plt.gca().set_aspect('equal', adjustable='box')
#contour = plt.contour(Xtop,Ytop,ztop,8,colors='k')
#plt.clabel(contour, fontsize=10, colors='k')
    plt.xticks()
    plt.yticks()
    nfont=15
    
    
    #plt.xlim(-adir0[0, 0], adir0[0, 0]*1)
    #plt.ylim(-adir0[0, 0] * 3**0.5/2*1.2, adir0[0, 0] * 3**0.5/2*1.2)
    #plt.xlim([-adir0[0, 0], -adir0[0, 0] + 3 * adir0[0, 0]])
    #plt.ylim([-adir0[0, 0]* 3**0.5, -adir0[0, 0] * 3**0.5 + 1.5 * 3**0.5 * adir0[0, 0]])



    plt.xlim([-0.9, -1 + 2.4])
    plt.ylim([-1* 3**0.5*0.8, -1 * 3**0.5 +  1.3 * 3**0.5 * 1])



    adir = adir0 
    #plt.plot([0.0, adir[0,0]], [0.0, adir[0,1]], "--", color="k", linewidth=1)
    #plt.plot([0.0, 0.0], [0.0, adir[1,1]], "--", color="k", linewidth=1)
    #plt.plot([adir[0,0], adir[0, 0]], [adir[0, 1], adir[0, 1]+ adir[1,1]], "--", color="k", linewidth=1)
    #plt.plot([0.0, adir[0, 0]], [adir[1,1], adir[1,1]+adir[0, 1]], "--", color="k", linewidth=1)
    
    plt.plot([0.0-0.5, 1.0-0.5], [0.0-3**0.5/2, 0-3**0.5/2], "--", color="k", linewidth=1)
    plt.plot([0.0-0.5, 0.5-0.5], [0.0-3**0.5/2, 3**0.5/2-3**0.5/2], "--", color="k", linewidth=1)
    plt.plot([1.0-0.5, 1+0.5-0.5], [0-3**0.5/2, 3**0.5/2-3**0.5/2], "--", color="k", linewidth=1)
    plt.plot([0.5-0.5, 1.5-0.5], [3**0.5/2-3**0.5/2, 3**0.5/2-3**0.5/2], "--", color="k", linewidth=1)



    #orig = [ 0.5 * adir0[0, 0], 3**0.5/2 * adir0[0, 0]]
    #plt.plot([0.0-orig[0], adir0[0,0]-orig[0]], [0.0-orig[1], 0.0-orig[1]], ':', color='k', linewidth=1)
    #plt.plot([0.0-orig[0], adir0[1, 0]-orig[0]], [0.0-orig[1], adir0[1,1]-orig[1]], ':', color='k', linewidth=1)
    #plt.plot([adir0[0,0]-orig[0], adir0[0, 0]+ adir0[1,0]-orig[0]], [0.0-orig[1], adir0[1,1]-orig[1]], ':', color='k', linewidth=1)
    #plt.plot([adir0[1,0]-orig[0], adir0[0, 0]+ adir0[1,0]-orig[0]], [adir0[1,1]-orig[1], adir0[1,1]-orig[1]], ':', color='k', linewidth=1)
    

    #orig = [ 0.5 * adir0[0, 0], 3**0.5/2 * adir0[0, 0]]
    #plt.plot([0.0-orig[0], adir0[0,0]-orig[0]], [0.0-orig[1], 0.0-orig[1]], ':', color='w', linewidth=1)
    #plt.plot([0.0-orig[0], adir0[1, 0]-orig[0]], [0.0-orig[1], adir0[1,1]-orig[1]], ':', color='w', linewidth=1)
    #plt.plot([adir0[0,0]-orig[0], adir0[0, 0]+ adir0[1,0]-orig[0]], [0.0-orig[1], adir0[1,1]-orig[1]], ':', color='w', linewidth=1)
    #plt.plot([adir0[1,0]-orig[0], adir0[0, 0]+ adir0[1,0]-orig[0]], [adir0[1,1]-orig[1], adir0[1,1]-orig[1]], ':', color='w', linewidth=1)
    
    cb=plt.colorbar(cset,fraction=0.04, pad=0.08, orientation='horizontal')
    tick_locator = ticker.MaxNLocator(nbins=5)  # colorbar上的刻度值个数
    cb.locator = tick_locator
    cb.ax.yaxis.major.formatter._useMathText = True
    cb.ax.yaxis.major.formatter._useOffset = True
    #cb.ax.yaxis.set_major_formatter(OOMFormatter(-3, "%1.1f")) #OOMFormatter(-3, "%1.1f"))
    #cb.set_ticks([0.0, 1.0, 2.0, 3.0 ])
#cb.set_ticklabels([-1.5, -1, -0.5, 0.0, 0.5, 1, 1.5, 2.0, 2.5])
    cb.update_ticks()
    #plt.clim(-19, 19)
    cb.ax.set_title('$\Delta_{z}$(r) (meV)', size=nfont+2)
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
    plt.savefig('fig-paper-sky-.png',dpi=500, transparent=True)
    return 

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
    plt.savefig("plot-pot-diag"+sflag+".png", dpi=500, transparent=True)
    return 

grid = np.zeros(2, dtype=int)

grid[:] = [100, 100]

grid2 = np.zeros(2, dtype=int)
grid2[:] = [10, 10]

print("grid:", grid[0], grid[1])
Delta_t, Delta_b, Delta_T = read_pot(grid, grid2)

#min_val = min(np.reshape(pot_2d_top_sm, (grid[0]*grid[1], -1)))


plot_pot(grid, grid2, Delta_t, Delta_b, Delta_T)
#print("max:", max(np.reshape(pot_2d_top_sm-min_val, (grid[0]*grid[1], -1))))
#plot_pot_diag(adir0, pot_diag, 'top')


#f=open("pot_diag-top.dat", 'w')
#for i1 in range(0, grid[0]):
#    np.savetxt(f, np.mat(pot_diag[i1, :]), fmt="%20.14E %20.14E")
#f.close()
#min_val = min(np.reshape(pot_2d_bot_sm, (grid[0]*grid[1], -1)))
#plot_pot(adir0, grid, pot_2d_bot_sm-min_val, 'bot')

#print("max:", max(np.reshape(pot_2d_top_sm-pot_2d_bot_sm, (grid[0]*grid[1], -1))))
#print("min:", min(np.reshape(pot_2d_top_sm-pot_2d_bot_sm, (grid[0]*grid[1], -1))))

#plot_pot(adir0, grid, pot_2d_top_sm-pot_2d_bot_sm, 'top-bot')
