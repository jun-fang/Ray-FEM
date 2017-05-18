import numpy as np
import math

omega = np.array([314.159265358979,
471.238898038469,
628.318530717959,
942.477796076938,
1570.79632679490,
2356.19449019235,
3141.59265358979])

omega = omega/np.pi

err_NPW_4 = np.array([ 0.00145017580508982,
0.000814197073212150,
0.000578722488460479,
0.000372982064475219,
0.000223777788733937,
0.000152180307288643,
0.000113806282747749])


err_NPW_6 = np.array([ 0.000632524668118165,
0.000352240614940233,
0.000251410390616387,
0.000167663373601132,
0.000100997451113749,
6.69132403759567e-05,
5.02536762107745e-05])


err_NPW_8 = np.array([ 0.000350460197569831,
0.000194739591859630,
0.000144018468438206,
9.47447579277406e-05,
5.67219180239109e-05,
3.76416235198704e-05,
2.82002992152050e-05])



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

p4, = plt.loglog(omega[:len(err_NPW_4)], 0.6*err_NPW_4[0]/((omega[:len(err_NPW_4)]/(omega[0]))**1), 
	label=r'$\mathcal{O}(\omega^{-1})$', color='r', linewidth=2, linestyle= 'solid', markersize=8.0, zorder=2)

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
plt.xlim(90, 1100)
# plt.ylim(10**(-5), 10**(-2)/2)
plt.ylim(0.8*err_NPW_8[-1], 1.2*err_NPW_4[0])
plt.tight_layout(pad=0.5)

fig.savefig('ex3_Constant_Gradient_Velocity_ExRay.pdf')
plt.show()

plt.close('all')



