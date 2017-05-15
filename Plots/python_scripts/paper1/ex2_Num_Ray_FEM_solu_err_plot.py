# Plot Nunerical Ray-FEM Solution Error 


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

NRayFEM_solu_err1 = np.array([ 0.000263551282082452,
0.000254140256937518,
0.000101494076854674,
0.000124117503972135,
6.64798443047013e-05,
6.66127507918252e-05,
6.73550207225830e-05,
5.83311590181842e-05])


NRayFEM_solu_err2 = np.array([9.59168285391563e-05,
6.87601603333895e-05,
4.91275290610534e-05,
4.98225006405050e-05,
3.76452744600659e-05,
4.11255002113623e-05,
2.87953174477915e-05,
2.49374764198247e-05])


import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


p1, = plt.loglog(omega, NRayFEM_solu_err1, label=r'$\Vert u_{\mathbf{d}_{\widetilde{\omega}}} - u_{ex}\Vert_{L^2(\Omega)} $',
           color='b', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

p2, = plt.loglog(omega, NRayFEM_solu_err2, label=r'$\Vert u_{\mathbf{d}_{\omega}} - u_{ex}\Vert_{L^2(\Omega)}$',
           color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# plt.loglog(omega, ERayFEM_solu_err, label=r'$\Vert u_{\mathbf{d}_{ex}} - u_{ex}\Vert_{L^2(\Omega)}$',
#           color='r', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

p3, = plt.loglog(omega, NRayFEM_solu_err2[0]*1.02/(omega/(omega[0]))**0.5, label=r'$\mathcal{O}(\omega^{-1/2})$', 
           color='r', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

first_legend = plt.legend(handles=[p1, p2], loc=1, ncol=1, frameon=False, fontsize=22)
ax = plt.gca().add_artist(first_legend)
plt.legend(handles=[p3], loc=3, ncol=1, frameon=False, fontsize=22)

# plt.title('Numerical Ray-FEM',fontsize=20)

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Error', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.ylim(1.5*1e-5, .75*1e-3)
plt.tight_layout(pad=0.1)


fig.savefig('ex2_Num_Ray_FEM_solu_err_plot.pdf')
plt.show()

plt.close('all')


