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

ray_err1 = np.array([ 0.000189278099950591,
0.000157195022244089,
0.000124668153613167,
0.000107999522353882,
7.27900749535485e-05,
6.09210035627259e-05,
4.99141377789885e-05,
4.09631347989203e-05,
3.40825950959180e-05,
3.05009734099384e-05,
2.71400036417817e-05,
2.39507224476902e-05])


ray_err2 = np.array([  4.59946750271904e-05,
3.24000965513920e-05,
3.25577045548322e-05,
2.02505621467061e-05,
1.65708293384361e-05,
1.33338733488614e-05,
1.12437633267906e-05,
7.10482380540406e-06,
5.62001495169629e-06,
5.54170668650595e-06,
5.33968175527795e-06,
3.30567234462795e-06])

exray_norm = np.array([3.96106984480185,
3.95561014571193,
3.94878553580807,
3.94469077731188,
3.94091100519617,
3.93854865004105,
3.93650127707747,
3.93486337971239,
3.93377144863233,
3.93299149810422,
3.93240653534113,
3.93195156438202])

ray_err1 = ray_err1*exray_norm
ray_err2 = ray_err2*exray_norm


import matplotlib.pyplot as plt
# %matplotlib inline

golden = 1.61803398875
width = 6
height = width/golden



fig = plt.figure(figsize=(width, height))


plt.loglog(omega, ray_err1, label=r'$\Vert \theta (\mathbf{d}_{\widetilde{\omega}}) - \theta_{ex}\Vert_{L^2(\Omega)}$', 
           fillstyle = 'right', color='b', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, ray_err2, label=r'$\Vert \theta (\mathbf{d}_{\omega}) - \theta_{ex}\Vert_{L^2(\Omega)}$',
           color='m', linewidth=2, linestyle='--', marker='.', markersize=8.0, zorder=2)

plt.loglog(omega, ray_err2[0]/(omega/(omega[0])), label=r'$\mathcal{O}(\omega^{-1})$', color='g', linewidth=2, linestyle='solid', markersize=8.0, zorder=2)

# plt.loglog(N_x**2, N_x**2 / 4.0e4, label=r' ', color='white', linewidth=0.0)

plt.legend(loc=0, ncol=1, frameon=False, fontsize=16)

# plt.title('Normalized run-time for inner loop')

plt.xlabel(r'$\omega$', fontsize=18)
plt.ylabel(r'Error', fontsize=18)

plt.gca().tick_params(labelsize=14)

plt.autoscale(True, 'both', True)
plt.tight_layout(pad=0.1)

fig.savefig('ex1_ray_err_plot.pdf')
plt.show()

plt.close('all')


