import numpy as np
import math

omega = np.array([376.991118430775,
502.654824574367,
753.982236861550,
1005.30964914873,
1256.63706143592,
1570.79632679490,
1884.95559215388])

err_NPW_4 = np.array([ 0.00258930982124353,
0.00215155912258392,
0.00159809319661259,
0.00132432620681086,
0.00120835699929618,
0.000922924362604902])


err_NPW_6 = np.array([ 0.00146232132508687,
0.00130361456590141,
0.000879405066841832,
0.000717086709329996,
0.000646749254357502])


err_NPW_8 = np.array([ 0.000975879449384578,
0.000861420339184053,
0.000626648187115122,
0.000543591774964586])



import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


p1, = plt.loglog(omega[:len(err_NPW_4)], err_NPW_4, label=r'NPW = 4',
        	color='b', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

p2, = plt.loglog(omega[:len(err_NPW_6)], err_NPW_6, label=r'NPW = 6', 
			color='g', linewidth=2, linestyle= '--', marker='o', markersize=8.0, zorder=2) 

p3, = plt.loglog(omega[:len(err_NPW_8)], err_NPW_8, label=r'NPW = 8', 
			color='k', linewidth=2, linestyle= '--', marker='o', markersize=8.0, zorder=2)

p4, = plt.loglog(omega[:len(err_NPW_4)], 0.7*err_NPW_4[0]/((omega[:len(err_NPW_4)]/(omega[0]))**0.5), label=r'$\mathcal{O}(\omega^{-1/2})$', 
			color='r', linewidth=2, linestyle= 'solid', markersize=8.0, zorder=2)

# p1, = plt.loglog(omega, err_NPW_4, label=r'$\Vert u_{\mathbf{d}_{ex}} - u_{ex}\Vert_{L^2(\Omega)}$',
#            color='g', linewidth=2, linestyle='--', marker='o', markersize=8.0, zorder=2)

# p2, = plt.loglog(omega, err_NPW_6, label=r'$\mathcal{O}(\omega^{-1})$', color='r', linewidth=2, 
# linestyle= '--', markersize=8.0, zorder=2)  # linestyle= 'solid'

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)
first_legend = plt.legend(handles=[p1, p2, p3], loc=1, ncol=1, frameon=False, fontsize=15)
ax = plt.gca().add_artist(first_legend)
plt.legend(handles=[p4], loc=3, ncol=3, frameon=False, fontsize=18)

# plt.legend(loc=0, ncol=1, frameon=False, fontsize=26)

# plt.title('Exact Ray-FEM',fontsize=20)

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Rel $L^2$ Err', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.xlim(0.9*omega[0], 1.1*omega[-1])
plt.ylim(0.8*err_NPW_8[-1], 1.2*err_NPW_4[0])
plt.tight_layout(pad=0.5)

fig.savefig('ex3_Constant_Gradient_Velocity.pdf')
plt.show()

plt.close('all')



