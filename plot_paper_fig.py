import basemap_shaded as bs
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from custom_cmap import truncate_colormap

fig = plt.figure(figsize=(5,10))
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1],
                            width_ratios=[1, 1])

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

cmap = plt.get_cmap('terrain')
# cmap = plt.get_cmap('gist_earth')
newcmap = truncate_colormap(cmap,0.35,0.9)
newcmap.set_under('w')

bs.plot_obs_domain(ax=ax1, cmap=newcmap)
bs.plot_petaluma_gap(ax=ax2,cmap=newcmap)

plt.show(block=False)
