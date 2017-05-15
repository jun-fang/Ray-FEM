import numpy as np
import math

omega = np.array([125.663706143592,
150.796447372310,
201.061929829747,
251.327412287183,
326.725635973339,
402.123859659494,
502.654824574367,
628.318530717959,
753.982236861551,
879.645943005143,
1005.30964914873,
1130.97335529233])

NRayFEM_solu_err1 = np.array([ 9.09346707650119e-05,
7.30188440559808e-05,
5.56195075944363e-05,
4.01398650789899e-05,
2.94613874774725e-05,
2.42913214175416e-05,
1.89335930032249e-05,
1.58221393837071e-05,
1.31414667570940e-05,
1.12424440539432e-05,
9.84802158754574e-06,
8.75334024485210e-06])


NRayFEM_solu_err2 = np.array([ 6.57088019089599e-05,
4.58491377058688e-05,
4.05049681131386e-05,
3.06728032529348e-05,
2.20784199890883e-05,
1.95470949161170e-05,
1.58790421719717e-05,
1.23922970956897e-05,
1.00104424775907e-05,
8.97108870850093e-06,
7.82315776374285e-06,
6.78051120872111e-06])

ERayFEM_solu_err = np.array([6.18886282382786e-05,
5.17804696619368e-05,
3.88982313939246e-05,
3.11712871565166e-05,
2.40213809955053e-05,
1.95418517394008e-05,
1.56597004242518e-05,
1.25369059551770e-05,
1.04510257422425e-05,
8.96091605062924e-06,
7.84390728568119e-06,
6.97507638869924e-06])


Ex_solu_norm = np.array([0.479700650321934,
0.479036217069733,
0.478205671965820,
0.477707343058336,
0.477247345940659,
0.476959847164174,
0.476710681202540,
0.476511348197756,
0.476378459412448,
0.476283538795265,
0.476212348301704,
0.476156977899733])

NRayFEM_solu_err1 = NRayFEM_solu_err1*Ex_solu_norm
NRayFEM_solu_err2 = NRayFEM_solu_err2*Ex_solu_norm
ERayFEM_solu_err = ERayFEM_solu_err*Ex_solu_norm


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

p3, = plt.loglog(omega, 0.9*NRayFEM_solu_err2[0]/(omega/(omega[0])), label=r'$\mathcal{O}(\omega^{-1})$', 
           color='r', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

first_legend = plt.legend(handles=[p1, p2], loc=1, ncol=1, frameon=False, fontsize=23)
ax = plt.gca().add_artist(first_legend)
plt.legend(handles=[p3], loc=3, ncol=1, frameon=False, fontsize=22)

#plt.title('Numerical Ray-FEM',fontsize=20)

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Error', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.ylim(2.5*1e-6, 1e-4)
plt.tight_layout(pad=0.1)


fig.savefig('ex1_Num_Ray_FEM_solu_err_plot.pdf')
plt.show()

plt.close('all')



