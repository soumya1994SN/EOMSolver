# EOMSolver
`EOMSolver` is a python package for solving the equations of motion governing fast flavor evolution of a dense neutrino gas. Please see `1DBoxsolver.py` and `2DBoxsolver.py` for more details.

`1DBoxsolver.py` : Solves Eq. (3.2.1) in the thesis as a function of time. 

`2DBoxsolver.py` : Solves Eq. (4.3.1) in the thesis as a function of time. 

`EOMSolver` requires `Intelpython` and the recent `SciPy` and `NumPy` libraries. The code is written in polarization vector language. The inputs that need to be supplied are 
- The three components of the polarization vectors for all velocity modes and spatial locations at time t = 0. 
- Spatial discretizations
- Time discretization.
- Velocity discretizations. 


