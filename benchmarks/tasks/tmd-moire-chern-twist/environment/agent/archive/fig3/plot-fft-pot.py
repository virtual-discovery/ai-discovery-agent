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
def plot_pot_diag(line_grid, pot_diag_0, pot_diag_1, pot_diag_2, pot_diag_3):
    fig, ax1 = plt.subplots()
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params( width=1.5, direction='in', length=8, pad=8)  

    nfont = 33
    ax1.tick_params(axis='x', which='minor', direction='in', length=4, labelsize=nfont)    
    fig.set_size_inches(8, 6)
    ax1.xaxis.major.formatter._useMathText = True
    plt.rcParams['mathtext.fontset'] = 'stix'         

#for j in range(0, nbnd):
    
    ax1.plot(pot_diag_0[1:, 0], pot_diag_0[1:, 1]*1000,  linewidth=2, linestyle="dashed", marker="o" , markersize=9, c="red",  label="1.25$^{\circ}$")
    ax1.plot(pot_diag_1[1:, 0], pot_diag_1[1:, 1]*1000, linewidth=2, linestyle="dashed", marker="^", markersize=9,  c="orange", label="1.54$^{\circ}$")
    ax1.plot(pot_diag_2[1:, 0], pot_diag_2[1:, 1]*1000,linewidth=2, linestyle="dashed", marker="s",markersize=9,   c="green",  label="1.70$^{\circ}$")
    ax1.plot(pot_diag_3[1:, 0], pot_diag_3[1:, 1]*1000,linewidth=2, linestyle="dashed", marker="D", markersize=9,  c="blue", label="3.15$^{\circ}$")
    #ax1.plot(pot_diag_4[:, 0], pot_diag_4[:, 1]*1000,linestyle="solid", linewidth=3, color="#2B3990",  label="3.15$^{\circ}$")

    #ax1.scatter(pot_diag_0[1:, 0], pot_diag_0[1:, 1]*1000, marker="o" , markersize=9, c="red",  label="1.25$^{\circ}$")
    #ax1.scatter(pot_diag_1[1:, 0], pot_diag_1[1:, 1]*1000,  marker="^", markersize=9,  c="orange", label="1.47$^{\circ}$")
    #ax1.scatter(pot_diag_2[1:, 0], pot_diag_2[1:, 1]*1000, marker="s",markersize=9,   c="green",  label="1.70$^{\circ}$")
    #ax1.scatter(pot_diag_3[1:, 0], pot_diag_3[1:, 1]*1000,  marker="D", markersize=9,  c="blue", label="3.15$^{\circ}$")
    #ax1.plot(pot_diag_4[:, 0], pot_diag_4[:, 1]*1000,linestyle="solid", linewidth=3, color="#2B3990",  labe

    ax1.set_xlabel("g/g$_{1}$", fontsize=nfont+2, fontname="STIXGeneral")
    ax1.set_ylabel("$\\Delta$v$_{H}$(g) (meV)", fontsize=nfont+2, fontname="STIXGeneral")
    ax1.legend(frameon=False, fontsize=nfont,labelspacing = 0.0, loc='upper right', handlelength=1, handletextpad=0.2)


    plt.xlim([0.5, 7])
    #plt.ylim([-0.5, 25])
    #plt.xticks([])
    #plt.xlim(0, adir0[0, 0] *1.5)
    #plt.text(, -4.5, 'MM',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")
    #plt.text(0.0, -4.5, 'MM',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")
    #plt.text(0.0, -4.5, 'MM',horizontalalignment='center',  fontsize=nfont, fontname="STIXGeneral")

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
    plt.savefig("plot-pot-fft-degree.png", dpi=500)
    return




line_grid = np.zeros(4, dtype=int)
line_grid = [10000, 10000, 10000,  10000 ]
nline = 10000

iflag = 0 
pot_diag_0 = np.zeros((nline, 2), dtype=float)
f1=open("1.25/pot_top_bot_G.dat", "r")
for i1 in range(0, nline):
    line=f1.readline() 
    pot_diag_0[i1, :]=[float(line.split()[i]) for i in range(0,2)]
f1.close() 
dx = pot_diag_0[1, 0] - pot_diag_0[0, 0]
#pot_diag_0[nline, 0] = pot_diag_0[nline-1, 0] + dx 
#pot_diag_0[nline, 1] = pot_diag_0[0, 1]
pot_diag_0[0,1] = 0


iflag = 1
pot_diag_1 = np.zeros((nline, 2), dtype=float)
f1=open("1.54/pot_top_bot_G.dat", "r")
for i1 in range(0, nline):
    line=f1.readline() 
    pot_diag_1[i1, :]=[float(line.split()[i]) for i in range(0,2)]
f1.close() 
dx = pot_diag_1[1, 0] - pot_diag_1[0, 0]
#pot_diag_1[nline, 0] = pot_diag_1[nline-1, 0] + dx 
#pot_diag_1[nline, 1] = pot_diag_1[0, 1]
pot_diag_1[0, 1] = 0



iflag = 2
f1=open("1.70/pot_top_bot_G.dat", "r")
pot_diag_2 = np.zeros((nline, 2), dtype=float)
for i1 in range(0, nline):
    line=f1.readline() 
    pot_diag_2[i1, :]=[float(line.split()[i]) for i in range(0,2)]
f1.close() 
dx = pot_diag_2[1, 0] - pot_diag_2[0, 0]
#pot_diag_2[nline, 0] = pot_diag_2[nline-1, 0] + dx 
#pot_diag_2[nline, 1] = pot_diag_2[0, 1]
pot_diag_2[0, 1] = 0 



iflag = 3
f1=open("3.15/pot_top_bot_G.dat", "r")
pot_diag_3 = np.zeros((nline, 2), dtype=float)
for i1 in range(0, nline):
    line=f1.readline() 
    pot_diag_3[i1, :]=[float(line.split()[i]) for i in range(0,2)]
f1.close() 
dx = pot_diag_3[1, 0] - pot_diag_3[0, 0]
#pot_diag_4[nline, 0] = pot_diag_4[nline-1, 0] + dx 
#pot_diag_4[nline, 1] = pot_diag_4[0, 1]
pot_diag_3[0, 1] = 0 

plot_pot_diag(line_grid, pot_diag_0, pot_diag_1, pot_diag_2,  pot_diag_3)
