import numpy as np
import math

omega = np.array([31.4159265358979,
62.8318530717959,
94.2477796076938,
125.663706143592,
188.495559215388,
251.327412287183,
314.159265358979])

NRayFEM_solu_err1 = np.array([ 0.000572648891767581,
0.000722865928822234,
0.000864110220427005,
0.000419176193721787,
0.000218112092091791,
0.000172378810810646,
0.000186612852486354])


NRayFEM_solu_err2 = np.array([ 0.000254791766422364,
0.000168749640873052,
0.000144488606792547,
0.000133444478123013,
0.000122936858862264,
0.000108502939013696,
9.24760276897857e-05])

ERayFEM_solu_err = np.array([ 1.93987862131576e-05,
9.87568976717513e-06,
6.55130518677237e-06,
4.93930497821691e-06,
3.29151693456906e-06,
2.46730294271270e-06,
1.97327245089172e-06])

Ex_solu_norm = np.array([2.20640163896561,
2.16646448922153,
2.15318377187808,
2.14654935800579,
2.13991890948041,
2.13660517271996,
2.13461740677809])

NRayFEM_solu_err1 = NRayFEM_solu_err1*Ex_solu_norm
NRayFEM_solu_err2 = NRayFEM_solu_err2*Ex_solu_norm
ERayFEM_solu_err = ERayFEM_solu_err*Ex_solu_norm



import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


plt.loglog(omega, NRayFEM_solu_err1, label=r'$\Vert u_{\widetilde{\omega},h}^{NR} - u_{ex}\Vert_{L^2} /{\Vert u_{ex}\Vert_{L^2}}$',
           color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, NRayFEM_solu_err2, label=r'$\Vert u_{\omega,h}^{NR} - u_{ex}\Vert_{L^2}/{\Vert u_{ex}\Vert_{L^2}}$',
           color='m', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, ERayFEM_solu_err, label=r'$\Vert u_{\omega,h}^{ER} - u_{ex}\Vert_{L^2}/{\Vert u_{ex}\Vert_{L^2}}$',
           color='r', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, 1/(omega/math.pi)/550, label=r'$\mathcal{O}(\omega^{-1})$', color='g', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=0, ncol=1, frameon=False, fontsize=14)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Ray-FEM Solution Error', fontsize=12)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.1)

fig.savefig('ex2_ray_solu_err_plot.pdf')
plt.show()

plt.close('all')


