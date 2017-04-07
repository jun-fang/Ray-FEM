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

timing_LowFreq = np.array([ 18.684,
31.671,
57.522,
72.086,
107.04,
139.11,
177.29,
213.87,
287.15,
315.41,
646.75,
490.08,
599.22,
776.20,
864.50,
943.94,
1254.5,
1392.8,
1383.9])


timing_HighFreq = np.array([ 7.767797,
11.69330,
16.92632,
22.92285,
42.91490,
41.97009,
57.93778,
64.64723,
81.93476,
98.76194,
127.5492,
150.4763,
171.8988,
181.0325,
232.2940,
302.6154,
289.2646,
341.9515,
389.1464])

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

plt.loglog(xlabels, timing_LowFreq, label='Low Frequency', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


# plt.plt.loglog(xlabels, timing_LowFreq, label='LowFrequency', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)
plt.loglog(size**2, timing_HighFreq, label='High Frequency', color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# plt.loglog(size**2, timing_gauss, label='Gaussian bumps', color='g', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(xlabels, (xlabels*np.log(xlabels)**2/(xlabels[0]*np.log(xlabels[0])**2))*timing_HighFreq[0]*0.95, label=r'$\mathcal{O}(N \log^2{N})$', color='k', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)
# #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# # plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=2, ncol=1, frameon=False, fontsize=14.85)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$N=n^2$', fontsize=18)
plt.ylabel('Nit', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.1)

fig.savefig('Low_freq_scaling.pdf')

plt.close('all')

# plt.show()


