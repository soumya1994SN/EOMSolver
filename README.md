# EOMSolver
`EOMSolver` is a python package for solving the equations of motion governing fast flavor evolution of a dense neutrino gas. Please see `1DBoxsolver.py` and `2DBoxsolver.py` for more details.

`1DBoxsolver.py` : Solves Eq. (3.2.1) in the thesis as a function of time. 

`2DBoxsolver.py` : Solves Eq. (4.3.1) in the thesis as a function of time. 

`EOMSolver` requires `Intelpython` and the recent `SciPy` and `NumPy` libraries. Especially, the `ode` and `fftpack` solvers under the `SciPy` library needs to be installed. The code is written in polarization vector language. The inputs that need to be supplied are :
- The initial condition for all the three components of the polarization vectors of all velocity modes at each spatial location. 
- Number of discretized bins in space ( ${N}_{sp}$)
- Number of discretized bins in time ( ${N}_{t}$)
- Number of discretized bins in velocity ( ${N}_{vel}$)

The output data represents $\mathsf{S}_{\vec{v}}[\vec{r}, t]$ at each $t, \vec{r}$ and $\vec{v}$ points. 

The data is generated as .npy files with shape $(N_{t}, N_{sp}, N_{vel}, 3)$ and $(N_{t}, N_{sp}, N_{sp}, N_{vel}, 3)$ for `1DBoxsolver.py` and `2DBoxsolver.py` respectively. In these .npy files the data is organized in an increasing order of each of the index. The indexes with shape $N_{t}, N_{sp}$ and $N_{sp}$ respectively indicate the data at each time division, spatial location and angular bin. The index of shape 3 represents the three components of the polarization vector in flavor space. For example, this index with value $0$, $1$ and $2$ respectively indicate the $\hat{\mathsf{e}}_1$, $\hat{\mathsf{e}}_2$ and $\hat{\mathsf{e}}_3$ components of the polarization vector for a specific space-time location and angular bin. 
