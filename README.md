# EOMSolver
`EOMSolver` is a python package for solving the equations of motion governing fast flavor evolution of a dense neutrino gas. Please see `1DBoxsolver.py` and `2DBoxsolver.py` for more details.

`1DBoxsolver.py` : Solves Eq. (3.2.1) in the thesis as a function of time. 

`2DBoxsolver.py` : Solves Eq. (4.3.1) in the thesis as a function of time. 

`EOMSolver` requires `Intelpython` and the recent `SciPy` and `NumPy` libraries. The code is written in polarization vector language. The inputs that need to be supplied are 
- The initial condition for all the three components of the polarization vectors of all velocity modes at each spatial location. 
- Number of discretized bins in space <img src="https://latex.codecogs.com/svg.image?\left(N_{sp}\right)" title="https://latex.codecogs.com/svg.image?\left(N_{sp}\right)" />
- Number of discretized bins in time ($N_t$)
- Number of discretized bins in velocity ($N_{vel}$)



