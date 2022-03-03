This directory contains the BFGS-Wolfe algorithm, the Standard BFGS-Armijo algorithm and the Standard BFGS-Wolfe algorithm for solving multiobjective optimization problems, described in

L. F. Prudente and D. R. Souza, A quasi-Newton method with Wolfe line searches for multiobjective optimization, technical report, 2022.

Date: March of 2022
License: This code is released under the GNU General Public License

This folder also contains the third-party free codes: 
1) software Algencan 3.1.1;
    -  E. G. Birgin and J. M. Martı́nez, Practical augmented Lagrangian methods for constrained optimization, SIAM, 2014.
    - https://www.ime.usp.br/~egbirgin/tango/

2) subroutines dcsrch and dcstep of Moré and Thuente.
    - J. J. Moré and D. J. Thuente, Line Search Algorithms with Guaranteed 
      Sufficient Decrease, ACM Trans. Math. Softw., 20 (1994), pp. 286–307.
    - http://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/

File main.f90 contains the main program. Modify myproblem.f90 routine to solve your own problem. Alternatively, set a test problem in main.f90 routine - see myproblem.f90.

Instructions:
—————————————

1) Go to folder and type 

make

2) Run typing

./MOPsolver

and see the output in the screen.
