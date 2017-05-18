import numpy as np
import math

omega = np.array([376.991118430775,
502.654824574367,
753.982236861550,
1005.30964914873,
1507.96447372310,
2010.61929829747,
3015.92894744620])

omega = omega/np.pi

err_NPW_4 = np.array([0.00270157770281540,
0.00219734201275802,
0.00170110076048316,
0.00141722055515402,
0.00116206600585611,
0.000834283365116496,
0.000760021129535664])


err_NPW_6 = np.array([ 0.00155179897461391,
0.00117260261040538,
0.000956857981728328,
0.000729819228791556,
0.000586008771195049,
0.000542605986265505,
0.000427214395423089])


err_NPW_8 = np.array([ 0.00103175668543789,
0.000777913307163148,
0.000617847788595371,
0.000492830804433761,
0.000408654829864075,
0.000354809212835210,
0.000264580813005265])



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

p4, = plt.loglog(omega[:len(err_NPW_4)], 0.75*err_NPW_4[0]/((omega[:len(err_NPW_4)]/(omega[0]))**0.65), 
	label=r'$\mathcal{O}(\omega^{-0.65})$', color='r', linewidth=2, linestyle= 'solid', markersize=8.0, zorder=2)

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

plt.xlabel(r'$\omega/\pi$', fontsize=18)
plt.ylabel(r'Rel $L^2$ Err', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
# plt.xlim(0.9*omega[0], 1.1*omega[len(err_NPW_4)-1])
plt.xlim(100, 1100)
plt.ylim(0.8*err_NPW_8[-1], 1.4*err_NPW_4[0])
plt.tight_layout(pad=0.5)

fig.savefig('ex1_NumRay_ConvRate.pdf')
plt.show()

plt.close('all')



