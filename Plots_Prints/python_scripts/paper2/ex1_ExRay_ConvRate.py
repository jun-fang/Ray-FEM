import numpy as np
import math

omega = np.array([376.991118430775,
502.654824574367,
753.982236861550,
1005.30964914873,
1507.96447372310,
2010.61929829747])

err_NPW_4 = np.array([ 0.00127617602407004,
0.000940066754055120,
0.000621635146193910,
0.000465273438145644,
0.000309519621286561,
0.000232549439626890])


err_NPW_6 = np.array([ 0.000568787774203152,
0.000418783338655105,
0.000276391902915176,
0.000206761799337624,
0.000137644360636887,
0.000103213473637043])


err_NPW_8 = np.array([ 0.000319120206092677,
0.000235864979624476,
0.000155349101652515,
0.000116330234214263,
7.74022989424810e-05,
5.80423659077803e-05])



import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


p1, = plt.loglog(omega, err_NPW_4, label=r'NPW = 4',
        	color='b', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

p2, = plt.loglog(omega, err_NPW_6, label=r'NPW = 6', 
			color='g', linewidth=2, linestyle= '--', marker='o', markersize=8.0, zorder=2) 

p3, = plt.loglog(omega, err_NPW_8, label=r'NPW = 8', 
			color='k', linewidth=2, linestyle= '--', marker='o', markersize=8.0, zorder=2)

p4, = plt.loglog(omega, 0.7*err_NPW_4[0]/(omega/(omega[0])), label=r'$\mathcal{O}(\omega^{-1})$', 
			color='r', linewidth=2, linestyle= 'solid', markersize=8.0, zorder=2)

# p1, = plt.loglog(omega, err_NPW_4, label=r'$\Vert u_{\mathbf{d}_{ex}} - u_{ex}\Vert_{L^2(\Omega)}$',
#            color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# p2, = plt.loglog(omega, err_NPW_6, label=r'$\mathcal{O}(\omega^{-1})$', color='r', linewidth=2, 
# linestyle= '--', markersize=8.0, zorder=2)  # linestyle= 'solid'

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)
first_legend = plt.legend(handles=[p1, p2, p3], loc=1, ncol=1, frameon=False, fontsize=14)
ax = plt.gca().add_artist(first_legend)
plt.legend(handles=[p4], loc=3, ncol=3, frameon=False, fontsize=20)

# plt.legend(loc=0, ncol=1, frameon=False, fontsize=26)

# plt.title('Exact Ray-FEM',fontsize=20)

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Error', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.xlim(0.9*omega[0], 1.1*omega[-1])
plt.ylim(0.8*err_NPW_8[-1], 1.2*err_NPW_4[0])
plt.tight_layout(pad=0.5)

fig.savefig('ex1_ExRay_ConvRate.pdf')
plt.show()

plt.close('all')



