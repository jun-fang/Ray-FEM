# Adaptive ray-FEM for solving high-frequency Helmholtz equation


We present a ray-based finite element method (ray-FEM) for the high-frequency Helmholtz equation in smooth media, whose basis are learned adaptively from the medium and source. 

The method requires a fixed number of grid points per wavelength to represent the wave field; moreover, it achieves an asymptotic convergence rate of  ![eq1](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BO%7D%28%5Comega%5E%7B-%5Cfrac%7B1%7D%7B2%7D%7D%29), where ![eq2](https://latex.codecogs.com/gif.latex?%5Comega) is the frequency parameter in the Helmholtz equation. The local basis are motivated by the geometric optics ansatz and are composed of polynomials modulated by plane waves propagating in a few dominant ray directions. The ray directions are learned by processing a low-frequency wave field that probes the medium with the same source. Once the local ray directions are extracted, they are incorporated into the local basis to solve the high-frequency Helmholtz equation. This process can be continued to further improve the approximations for both local ray directions and high-frequency wave fields iteratively. 

Finally, a fast solver is developed for solving the resulting linear system with an empirical complexity ![eq3](https://latex.codecogs.com/gif.latex?%5Cmathcal%7BO%7D%28%5Comega%5Ed%29) up to a poly-logarithmic factor. Numerical examples in 2D are presented to corroborate the claims.

## Papers

[Learning dominant wave directions for plane wave methods for high-frequency Helmholtz equations](https://link.springer.com/content/pdf/10.1186%2Fs40687-017-0098-9.pdf)

[A hybrid approach to solve the high-frequency Helmholtz equation with source singularity](https://arxiv.org/pdf/1710.02307.pdf)
