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

err_NPW_4 = np.array([ 0.00132817027775784,
0.000960374998031966,
0.000635954223250241,
0.000475999785279841,
0.000316521866452413,
0.000237678129843111,
0.000158657388654266])


err_NPW_6 = np.array([ 0.000585886057149803,
0.000429387568819536,
0.000282937622168569,
0.000211563030284355,
0.000140782318030535,
0.000105538174280112,
7.03701258054529e-05])


err_NPW_8 = np.array([ 0.000327867785064894,
0.000241418322140888,
0.000159022720684501,
0.000119036836395143,
7.91748238310790e-05,
5.93587070008570e-05,
3.95715164048433e-05])



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

plt.xlabel(r'$\omega/\pi$', fontsize=18)
plt.ylabel(r'Rel $L^2$ Err', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
# plt.xlim(0.9*omega[0], 1.1*omega[-1])
plt.xlim(100, 1100)
plt.ylim(0.8*err_NPW_8[-1], 1.2*err_NPW_4[0])
plt.tight_layout(pad=0.5)

fig.savefig('ex1_ExRay_ConvRate.pdf')
plt.show()

plt.close('all')



