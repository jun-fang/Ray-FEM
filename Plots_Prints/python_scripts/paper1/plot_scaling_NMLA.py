import numpy as np

size = np.array([45.0,
55.0,
65.0,
75.0,
85.0,
95.0,
105.0,
110.0,
120.0,
135.0,
145.0,
155.0,
160.0,
170.0,
185.0,
195.0,
205.0,
215.0,
225.0])

timing_NMLA = np.array([ 2.9972,
2.6158,
8.4915,
9.0110,
8.9726,
9.0912,
9.5124,
9.2398,
9.6367,
10.143,
9.2135,
9.9841,
20.423,
22.343,
22.594,
25.456,
24.850,
21.250,
21.9610])


# timing_LowFact = np.array([12.69,
# 21.37,
# 29.35,
# 41.82,
# 51.83,
# 69.76,
# 81.53,
# 104.58,
# 112.07,
# 124.82,
# 178.83,
# 211.56,
# 216.16,
# 289.90,
# 324.63,
# 282.74,
# 417.06,
# 461.85,
# 469.43 ])

# timing_gauss = np.array([ 1.218088277111111,
# 5.649571164,
# 8.722732030944444,
# 20.264415925555554,
# 34.40783213327778,
# 56.12274424983333,
# 70.80124191294445,
# 94.6494892356111,
# 148.3274009705])


import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))

xlabels = size**2;

plt.loglog(xlabels, timing_NMLA, label='Probing', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


# plt.plt.loglog(xlabels, timing_LowFreq, label='LowFrequency', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
#plt.loglog(size**2, timing_LowFact, label='Setup', color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# plt.loglog(size**2, timing_gauss, label='Gaussian bumps', color='g', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

# plt.loglog(xlabels, (xlabels*np.log(xlabels)**3/(xlabels[0]*np.log(xlabels[0])**3))*timing_LowFreq[0]*0.95, label=r'$\mathcal{O}(N \log^2{N})$', color='k', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)
# #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.loglog(xlabels, (xlabels*np.log(xlabels)/(xlabels[0]*np.log(xlabels[0])))*timing_NMLA[0]*0.95, label=r'$\mathcal{O}(N \log{N})$', color='r', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# # plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=2, ncol=1, frameon=False, fontsize=14.85)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$N=n^2$', fontsize=18)
plt.ylabel('Time [s]', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.2)

fig.savefig('Probing_scaling.pdf')

plt.close('all')

# plt.show()


