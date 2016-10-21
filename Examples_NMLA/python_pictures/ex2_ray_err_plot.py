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

ray_err1 = np.array([ 0.0483824686873896,
0.0286653601604059,
0.0278732156560360,
0.0346900764636998,
0.0285370726600453,
0.0161890874162595,
0.0155720147588110,
0.0135635762467273])


ray_err2 = np.array([   0.0166582871216415,
0.0114564464657079,
0.00844460800684360,
0.00725494435276593,
0.00582847267003839,
0.00495885945203899,
0.00441989900622282,
0.00384983177207796])




import matplotlib.pyplot as plt

golden = 1.61803398875
width = 6
height = width/golden

fig = plt.figure(figsize=(width, height))


plt.loglog(omega, ray_err1, label='$\Vert \\theta_{\widetilde{\omega},h} - \\theta_{ex}\Vert_{L^2(\Omega)}$', 
           fillstyle = 'right', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, ray_err2, label='$\Vert \\theta_{\omega,h} - \\theta_{ex}\Vert_{L^2(\Omega)} $',
           color='m', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, ray_err2[0]/(omega/(omega[0]))**0.5, label=r'$\mathcal{O}(\omega^{-1/2})$', color='g', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=0, ncol=1, frameon=False, fontsize=16)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Error', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.1)

fig.savefig('ex2_ray_err_plot.pdf')
plt.show()

plt.close('all')



