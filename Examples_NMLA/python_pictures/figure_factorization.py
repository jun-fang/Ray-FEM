import numpy as np

size = np.array([200.0,
        400.0,
        500.0,
        800.0,
        1000.0,
        1250.0,
        1400.0,
        1600.0,
        2000.0])

timing_factorization = np.array([ 0.77352955,
2.966673073,
4.260938869,
12.292629506,
17.698088652,
28.585287729,
35.76267899,
52.695841569,
79.456514489])


timing_build = np.array([   2.132851,
  6.304390,
  9.777168,
 24.985038,
 40.015493,
 61.975112,
 77.392156,
 96.510941,
150.347302])


import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


plt.loglog(size**2, timing_factorization, label='Factorization', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(size**2, timing_build, label='Assembly', color='m', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(size**2, (size**2)/(size[0]**2)*timing_factorization[0]*0.8, label=r'$\mathcal{O}(N)$', color='g', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=2, ncol=1, frameon=False, fontsize=14.85)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$N=n^2$', fontsize=18)
plt.ylabel(r'$t[s]$', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.1)

fig.savefig('scaling_factorization.pdf')

plt.close('all')

# plt.show()


