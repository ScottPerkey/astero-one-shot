#####THIS IS WHERE PANEL CODE BEGINS
from pudb import set_trace as pause
import pandas
import pylab as plt
import pickle
import numpy as np
import matplotlib. pyplot as plt
infile = open('data.pickle','rb')
pc = pickle.load(infile)
import matplotlib.gridspec as GridSpec
pause()

dfs = pc[0]
for p in pc[1:]:
    dfs = dfs.append(p)

# JCZ 110322
# have to make feh to be numerical
# '-1e-1'
dfs['feh'] = dfs['feh'].map(lambda x : x.replace('d', 'e'))
dfs['feh'] = dfs['feh'].astype(float)

first = dfs.loc[dfs['alpha'] == 2.0]
second = dfs.loc[dfs['alpha'] == 1.7]



panel_1=first.loc[first['mass'] == 0.8]
panel_2=first.loc[first['mass'] == 1.0]
panel_3=first.loc[first['mass'] == 1.2]
panel_4=second.loc[second['mass'] == 0.8]
panel_5=second.loc[second['mass'] == 1.0]
panel_6=second.loc[second['mass'] == 1.2]

# JCZ 110322
# only need to use either subplots OR gridspec. so commented out subplots
# fig, _ = plt.subplots(nrows=3, ncols=2, sharex= True,sharey= True)
# JCZ 110322
# if doing this, then need to then define a figure. so added this:
fig = plt.figure()
gs = fig.add_gridspec(nrows=3, ncols=2)

ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])
ax3 = fig.add_subplot(gs[3])
ax4 = fig.add_subplot(gs[4])
ax5 = fig.add_subplot(gs[5])
   

# JCZ 110322
# this of course needs to have all of the metallicities, not just 0 and 0.2.
fh= [-0.4, 0.0, 0.2]

# JCZ 110322
# commented this out entirely
# for feh in fh:
#     label= '[Fe/H]= {:3.2f}'.format(feh)
#     # JCZ 110322
#     # added print statements
#     # print(panel_6.loc[panel_6['feh'] == feh]['small_dnu'])
#     # print(panel_6.loc[panel_6['feh'] == feh]['dnu_obs'])
#     # JCZ 110322
#     # don't plot the points less than 0
#     panel_6 = panel_6.loc[panel_6['small_dnu'] > 0]
#     ax5.scatter(panel_6.loc[panel_6['feh'] == feh]['small_dnu'], panel_6.loc[panel_6['feh']== feh]['dnu_obs'], label=label) # JCZ 110322 added label=label

    # JCZ 110322
    # this makes all the axis plots
for ax, panel in zip([ax0, ax1, ax2, ax3, ax4, ax5], [panel_1, panel_2, panel_3, panel_4, panel_5, panel_6]):
    for feh in fh:
        label= '[Fe/H]= {:3.2f}'.format(feh)
        # JCZ 110322
        # don't plot the points less than 0

        # pick out the points that correspond to the pre-main sequence... this will need to actually be chosen by going back to the models and finding out when hydrogen burning begins. but this is a first step
        # actually, the numax values don't seem to be in an order.
        # pause()
        panel = panel.loc[panel['small_dnu'] > 0]
        # print(panel.loc[panel['feh'] == feh]['small_dnu'])
        # pause()
        ax.scatter(panel.loc[panel['feh'] == feh]['small_dnu'], panel.loc[panel['feh']== feh]['dnu_obs'], label=label) # JCZ 110322 added label=label



plt.tight_layout()
# JCZ 110322
# adding legend
plt.legend()
# ax.autoscale()

plt.savefig('pma.png', format='png'); plt.legend()
