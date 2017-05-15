import numpy as np

size = 10*np.array([45.0,
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

# timing_HighFreq = np.array([ 7.767797,
# 11.69330,
# 16.92632,
# 22.92285,
# 42.91490,
# 41.97009,
# 57.93778,
# 64.64723,
# 81.93476,
# 98.76194,
# 127.5492,
# 150.4763,
# 171.8988,
# 181.0325,
# 232.2940,
# 302.6154,
# 289.2646,
# 341.9515,
# 389.1464])


timing_LowFreq = np.array([ 24.10 ,
39.68 ,
60.17 ,
73.05 ,
95.33 ,
149.93,
165.38,
195.02,
265.58,
296.84,
385.56,
474.34,
516.32,
659.86,
762.89,
750.04,
990.06 ,
1124.23,
1423.43  ]);


timing_LowFact = np.array([ 13.71 ,
21.07 ,
28.78 ,
41.41 ,
52.45 ,
69.81 ,
82.49 ,
97.01 ,
108.30,
125.06,
171.96,
204.89,
215.91,
285.84,
327.01,
273.86,
419.73,
480.92 ,
463.34  ])

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

plt.loglog(xlabels, timing_LowFreq, label='Solve', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


# plt.plt.loglog(xlabels, timing_LowFreq, label='LowFrequency', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
plt.loglog(size**2, timing_LowFact, label='Setup', color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# plt.loglog(size**2, timing_gauss, label='Gaussian bumps', color='g', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(xlabels, (xlabels*np.log(xlabels)**3/(xlabels[0]*np.log(xlabels[0])**3))*timing_LowFreq[0]*0.95, label=r'$\mathcal{O}(N \log^3{N})$', color='k', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)
# #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.loglog(xlabels, (xlabels*np.log(xlabels)/(xlabels[0]*np.log(xlabels[0])))*timing_LowFact[0]*1.1, label=r'$\mathcal{O}(N\log{N})$', color='r', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# # plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=2, ncol=1, frameon=False, fontsize=14.85)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$N=n^2$', fontsize=18)
plt.ylabel('Time [s]', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.2)

fig.savefig('Low_freq_scaling_heterogeneous.pdf')

plt.close('all')

# plt.show()


