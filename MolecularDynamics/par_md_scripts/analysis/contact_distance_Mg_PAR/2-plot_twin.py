import numpy as np
import matplotlib.pyplot as plt
import matplotlib

custom_style = {'axes.axisbelow': True,
 'axes.edgecolor': 'black',
 'axes.facecolor': '#EAEAF2',
 'axes.grid': False,
 'axes.labelcolor': '.15',
 'axes.linewidth': 0,
 'font.family': 'Arial',
 'grid.color': 'gray',
 'grid.linestyle': '--',
 'image.cmap': 'Greys',
 'legend.frameon': False,
 'legend.numpoints': 1,
 'legend.scatterpoints': 1,
 'lines.solid_capstyle': 'round',
 'pdf.fonttype': 42,
 'text.color': '.15',
 'xtick.color': '.15',
 'xtick.direction': 'in',
 'xtick.major.size': 0,
 'xtick.minor.size': 0,
 'ytick.color': '.15',
 'ytick.direction': 'in',
 'ytick.major.size': 0,
 'ytick.minor.size': 0}

#sns.set_style("dark", rc=custom_style)

import matplotlib
font = {'family' : 'arial',
  'weight' : 'medium',
  'size'  : 20}
matplotlib.rc('font', **font)
plt.get_cmap("jet")
matplotlib.rcParams['axes.linewidth'] = 2 #set the values globally                               
matplotlib.rcParams['xtick.major.size'] = 4
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['xtick.minor.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 2
matplotlib.rcParams['ytick.minor.width'] = 2

ts = 2e-6*5000*1e-3

c = np.loadtxt("data_min_dists.dat")[1:]
ee = np.loadtxt("../data/end2end_dist_eq1-LR4.dat")[1:]
rg = np.loadtxt("../data/rgyr_eq1-LR4.dat")
t = ts*np.arange(1,len(c)+1)

fig, ax = plt.subplots(figsize=(7,5.))

w = 1

ax.plot(t, ee, lw=w, color='lime', label=r'R$_{EE}$', zorder=1)
ax.plot(t, rg, lw=w, color='dodgerblue', label=r'R$_{g}$', zorder=2)

ax.set_ylabel(r"Length ($\mathrm{\AA}$)", )

ax2 = ax.twinx()
ax2.plot(t, c, lw=w, color='red', zorder=0)

ax2.set_ylabel(r'Mg$^{2+}$--PAR contact distance ($\mathrm{\AA}$)', color='r')
ax2.tick_params('y', colors='r')

ax.legend(loc='upper left', ncol=2, handlelength=1, labelspacing=0.25)
#ax2.legend(loc='upper right')

ax.tick_params(which='major',direction="out",length=6,right=False,top=False,labelsize=18, width=1.5)
ax2.tick_params(which='major',direction="out",length=6,right=True,top=False,labelsize=18, width=1.5)

plt.xlim(0)
plt.xlabel(r"Time ($\mathrm{\mu}$s)")

plt.savefig("compare_15mer__contact_rg_ee.png", dpi=300, bbox_inches='tight')
