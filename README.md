# Ray_FEM


We present a ray-based finite element method (ray-FEM) for the high-frequency Helmholtz equation in smooth media, whose basis are learned adaptively from the medium and source. 

The method requires a fixed number of grid points per wavelength to represent the wave field; moreover, it achieves an asymptotic convergence rate of $\mathcal{O}(\omega^{-\frac{1}{2}})$, where $\omega$ is the frequency parameter in the Helmholtz equation. The local basis are motivated by the geometric optics ansatz and are composed of polynomials modulated by plane waves propagating in a few dominant ray directions. The ray directions are learned by processing a low-frequency wave field that probes the medium with the same source. Once the local ray directions are extracted, they are incorporated into the local basis to solve the high-frequency Helmholtz equation. This process can be continued to further improve the approximations for both local ray directions and high-frequency wave fields iteratively. 

Finally, a fast solver is developed for solving the resulting linear system with an empirical complexity $\mathcal{O}(\omega^d)$ up to a poly-logarithmic factor. Numerical examples in 2D are presented to corroborate the claims.
