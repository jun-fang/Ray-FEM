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


timing_HighFreq = np.array([ 1.176772e+01 ,
1.761384e+01 ,
2.845541e+01 ,
3.497920e+01 ,
4.249018e+01 ,
5.997056e+01 ,
7.392939e+01 ,
9.390443e+01 ,
1.074052e+02 ,
1.263909e+02 ,
1.367733e+02 ,
1.868517e+02 ,
2.031273e+02 ,
2.381643e+02 ,
2.584227e+02 ,
3.247313e+02 ,
3.500849e+02 ,
4.242830e+02 ,
4.807868e+02  ]);


timing_HighFact = np.array([   6.12 ,
  8.54 ,
 11.57 ,
 17.72 ,
 21.72 ,
 26.46 ,
 32.60 ,
 37.03 ,
 41.53 ,
 51.51 ,
 61.32 ,
 70.73 ,
 82.07 ,
 102.43,
 112.95,
 119.90,
 136.93 ,
 150.46 ,
 168.30   ])

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

plt.loglog(xlabels, timing_HighFreq, label='Solve', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


# plt.plt.loglog(xlabels, timing_HighFreq, label='HighFrequency', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
plt.loglog(size**2, timing_HighFact, label='Setup', color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# plt.loglog(size**2, timing_gauss, label='Gaussian bumps', color='g', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(xlabels, (xlabels*np.log(xlabels)/(xlabels[0]*np.log(xlabels[0])))*timing_HighFreq[0]*1.05, label=r'$\mathcal{O}(N \log{N})$', color='k', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)
# #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.loglog(xlabels, (xlabels/(xlabels[0]))*timing_HighFact[0]*1.1, label=r'$\mathcal{O}(N)$', color='r', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# # plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=2, ncol=1, frameon=False, fontsize=14.85)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$N=n^2$', fontsize=18)
plt.ylabel('Time [s]', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.2)

fig.savefig('High_freq_scaling_heterogeneous.pdf')

plt.close('all')

# plt.show()


