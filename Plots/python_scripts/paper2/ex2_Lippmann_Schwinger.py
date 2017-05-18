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


err_NPW_4 = np.array([ 0.00216701147023616,
0.00176295146124669,
0.00130337855620765,
0.000919437663223020,
0.000673609558856852,
0.000595467258613511,
0.000440141889521627])


err_NPW_6 = np.array([ 0.00149085703643568,
0.00113989914906664,
0.000769409511702538,
0.000543214611743842,
0.000463120362868336,
0.000356876117029469,
0.000265312921615620])


err_NPW_8 = np.array([ 0.000993973425802366,
0.000727275400853260,
0.000505446079854768,
0.000405825911503581,
0.000302325437873150,
0.000238027167967451,
0.000198439997137332])



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

p4, = plt.loglog(omega, 0.8*err_NPW_4[0]/((omega/(omega[0]))**0.75), 
	label=r'$\mathcal{O}(\omega^{-0.75})$', color='r', linewidth=2, linestyle= 'solid', markersize=8.0, zorder=2)

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

plt.xlabel(r'$\omega/ \pi$', fontsize=18)
plt.ylabel(r'Rel $L^2$ Err', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.xlim(100, 1100)
plt.ylim(0.8*err_NPW_8[-1], 1.2*err_NPW_4[0])
plt.tight_layout(pad=0.5)

fig.savefig('ex2_Lippmann_Schwinger.pdf')
plt.show()

plt.close('all')



