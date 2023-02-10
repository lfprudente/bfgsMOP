This directory contains the BFGS-Wolfe algorithm, the Standard BFGS-Armijo algorithm and the Standard BFGS-Wolfe algorithm for solving multiobjective optimization problems described in the paper:

L. F. Prudente and D. R. Souza, [A quasi-Newton method with Wolfe line searches for multiobjective optimization](https://link.springer.com/article/10.1007/s10957-022-02072-5), *Journal of Optimization Theory and Applications* 194, pp. 1107-1140, 2022.
  

- MOPsolverBFGS.f90: routine containing the BFGS-Wolfe algorithm (a BFGS algorithm with the Hessian approximations updated at each iteration and step sizes satisfying the Wolfe condition)
- MOPsolverStBFGSArmijo.f90: routine containing the Standard BFGS-Armijo algorithm (a BFGS algorithm with a cautious update to the Hessian approximations and step sizes satisfying the Armijo condition)
- MOPsolverStBFGSWolfe.f90: routine containing the Standard BFGS-Wolfe algorithm (a BFGS algorithm with a cautious update to the Hessian approximations and step sizes satisfying the Wolfe condition)

This folder also contains the third-party free codes: 
- software Algencan 3.1.1
    -  E. G. Birgin and J. M. Martı́nez, Practical augmented Lagrangian methods for constrained optimization, SIAM, 2014.
    - https://www.ime.usp.br/~egbirgin/tango/
    - Algencan is used to compute the search directions; see innersolver.f90.

- subroutines dcsrch and dcstep of Moré and Thuente
    - J. J. Moré and D. J. Thuente, Line Search Algorithms with Guaranteed 
      Sufficient Decrease, ACM Trans. Math. Softw., 20 (1994), pp. 286–307.
    - http://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/
    - These subroutines are used as the inner solver of lsvecopt.f90 which computes a step size satisfying the (vector) Wolfe conditions.

Instructions:
-------------

File main.f90 contains the main program where you can choose the algorithm to be used.
Modify myproblem.f90 routine to solve your own problem. Alternatively, set a test problem in main.f90 routine; see myproblem.f90.

The codes are written in Fortran 90. Users need to install gfortran.

1) In the terminal, go to the folder and type:

    make

2) Run typing:

    ./MOPsolver

and see the output in the screen.

- out: outer iteration number
- |theta|: optimality measure 
- LS: flag of the line search routine to compute the step size (0 means success)
- IS: flag of the inner solver routine to compute the search direction (0 means success)
- #evalf: number of function evaluations
- #evalg: number of gradient evaluations
