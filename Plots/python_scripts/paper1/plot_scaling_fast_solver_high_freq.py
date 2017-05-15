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


timing_HighFreq = np.array([ 9.677481e+00 ,
1.197494e+01 ,
1.752305e+01 ,
2.346231e+01 ,
3.457710e+01 ,
4.967041e+01 ,
5.593979e+01 ,
7.619850e+01 ,
8.992520e+01 ,
1.001875e+02 ,
1.203322e+02 ,
1.388838e+02 ,
1.823083e+02 ,
1.953192e+02 ,
2.363119e+02 ,
2.921830e+02 ,
3.026785e+02 ,
3.399758e+02 ,
4.157844e+02]);



timing_HighFact = np.array([ 5.37,
 8.67,
12.83,
17.21,
21.91,
26.65,
30.97,
37.12,
43.63,
52.93,
64.71,
75.19,
84.40,
98.25,
109.15,
122.02,
128.57,
148.78,
158.01 ])

import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))

xlabels = size**2;

plt.loglog(xlabels, timing_HighFreq, label='Solve', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


# plt.plt.loglog(xlabels, timing_LowFreq, label='LowFrequency', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
plt.loglog(size**2, timing_HighFact, label='Setup', color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# plt.loglog(size**2, timing_gauss, label='Gaussian bumps', color='g', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(xlabels, (xlabels*np.log(xlabels)/(xlabels[0]*np.log(xlabels[0])))*timing_HighFreq[0]*0.95, label=r'$\mathcal{O}(N \log{N})$', color='k', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)
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

fig.savefig('high_freq_scaling_homogeneous.pdf')

plt.close('all')

# plt.show()


