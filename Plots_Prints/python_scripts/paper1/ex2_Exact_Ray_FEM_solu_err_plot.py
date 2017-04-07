import numpy as np
import math

omega = np.array([31.4159265358979,
62.8318530717959,
94.2477796076938,
125.663706143592,
188.495559215388,
251.327412287183,
376.991118430775,
502.654824574367])


ERayFEM_solu_err = np.array([9.91759008975895e-06,
4.55860969355428e-06,
2.97401903433026e-06,
2.16959724330724e-06,
1.42234451363623e-06,
1.05819207683400e-06,
7.00059070412104e-07,
5.22658917351204e-07])


import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


p1, = plt.loglog(omega, ERayFEM_solu_err, label=r'$\Vert u_{\mathbf{d}_{ex}} - u_{ex}\Vert_{L^2(\Omega)}$',
           color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

p2, = plt.loglog(omega, ERayFEM_solu_err[0]*0.9/(omega/(omega[0])), label=r'$\mathcal{O}(\omega^{-1})$', color='r', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)


first_legend = plt.legend(handles=[p1], loc=1, ncol=1, frameon=False, fontsize=25)
ax = plt.gca().add_artist(first_legend)
plt.legend(handles=[p2], loc=3, ncol=1, frameon=False, fontsize=22)


# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

# plt.legend(loc=0, ncol=1, frameon=False, fontsize=20)

# plt.title('Exact Ray-FEM',fontsize=20)

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Error', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.1)

fig.savefig('ex2_Exact_Ray_FEM_solu_err_plot.pdf')
plt.show()

plt.close('all')



