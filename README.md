# Ray_FEM


We present a ray-based finite element method (ray-FEM) by learning basis adaptive to the underlying high-frequency Helmholtz equation in smooth media. 

Based on the geometric optics ansatz of the wave field, we learn local dominant ray directions by probing the medium using low-frequency waves with the same source. Once local ray directions are extracted, they are incorporated into the finite element basis to solve the high-frequency Helmholtz equation. This process can be continued to further improve approximations for both local ray directions and the high-frequency wave field iteratively.

The method requires a fixed number of grid points per wavelength to represent the wave field and achieves an asymptotic convergence rate of $\omega^{-1/2}$ as the frequency $\omega \rightarrow \infty$. 

A fast solver is developed for the resulting linear system with an empirical complexity $O(\omega^d)$ up to a poly-logarithmic factor. Numerical examples in 2D are presented to corroborate the claims.

