import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Getting ts for time
tmp = np.loadtxt('data_10/monomer1.dat')
l = len(tmp)
t_total = 5000 * 2e-6 * l
frags = 10000
ts = 1. * t_total / frags
print(ts,t_total)
t_yaxis = np.arange(1000) * ts

data_monomers_mean = []
data_monomers_std = []
for m in range(1,16):
    x = np.loadtxt("data_10/monomer"+str(m)+".dat")
    x_split = np.array_split(x, frags)
    x_split_mean = np.array([np.mean(part) for part in x_split])
    x_split_std = np.array([np.std(part) for part in x_split])
    data_monomers_mean.append(x_split_mean)
    data_monomers_std.append(x_split_std)

data_monomers_mean = np.array(data_monomers_mean)
data_monomers_std = np.array(data_monomers_std)
print(data_monomers_mean.shape)
print(data_monomers_std.shape)

y_min = np.floor(np.min(data_monomers_mean-data_monomers_std))
y_max = np.ceil(np.max(data_monomers_mean+data_monomers_std))
print(y_min,y_max)

# Set thicker line width and tick width globally
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4

#fig, axes = plt.subplots(frags, 1, sharex=True, figsize=(5, 10))
plt.figure(figsize=(4,10))
ax = plt.gca()

cmap_style = 'magma'
imshow_obj = plt.imshow(data_monomers_mean.T, cmap=cmap_style, aspect='auto', origin='lower', 
		vmin=y_min, vmax=y_max, extent=(1,22,0,t_total))

plt.xlabel(r'PAR residue index', fontsize=20)
plt.ylabel(r'Time (ns)', fontsize=20, labelpad=5)

cbar = plt.colorbar(imshow_obj, label='Monomer contacts', aspect=35)
cbar.set_label(r'Number of PAR residues within 10 $\mathrm{\AA}$', fontsize=20, labelpad=9)  # Increase colorbar label font size
cbar.ax.tick_params(size=6,width=2)
cbar.outline.set_linewidth(2)  # Increase the linewidth of colorbar boundary
cbar.ax.tick_params(axis='y', labelsize=18)  # Increase colorbar tick label size
colorbar_min = y_min
colorbar_max = y_max 
# Set the colorbar limits directly on the imshow object
imshow_obj.set_clim(vmin=colorbar_min, vmax=colorbar_max)
cbar.mappable.set_clim(colorbar_min, colorbar_max)

plt.gca().spines['top'].set_linewidth(2)
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['right'].set_linewidth(2)

plt.xticks([1, 8, 16, 22], fontsize=18)
yticks = np.floor(np.linspace(0,t_total,6))
plt.yticks(yticks, fontsize=18)

plt.gca().tick_params(width=2)
plt.gca().tick_params(axis='x', length=6, width=2,labelsize=16)
plt.gca().tick_params(axis='y', length=6, width=2,labelsize=16)
plt.xticks([1,8,16,22])
#plt.yticks([])

#plt.gca().set_aspect('equal')  # Equal aspect ratio for a circular plot
plt.gca().set_aspect('auto')  # Adjust the aspect ratio

plt.savefig('15mer_100mM_NaCl_cluster10A_heatmap.pdf', dpi=300,bbox_inches = "tight")
