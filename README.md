# Mixed-Integer-ALADIN
This code uses the open-source tool CasADi:  
https://web.casadi.org/

Input problem must be partitioned:  
**obj_funs:** cell array of CasADi functions  
**cons_funs:** cell array of CasADi functions  
**A_mats:** cell array of matrices  
**b:** vector  
**x0:** cell array of vectors  
**lam0:** vector  
**lbx:** cell array of vectors  
**ubx:** cell array of vectors  
**discrete:** cell array of vectors  
**params:** struct  
**opts:** struct  

Params will assume default values if left empty. Possible params to set:  
**Sig:** cell array of matrices. Used to scale variables. Default: identity matrix  
**eps:** termination tolerance. Default: 10^-3  
**maxiter:** iteration limit. Default: 300  
**rho:** algorithm parameter analogous to step size. Default: 0.1obj(x0).  
**mu:** algorithm parameter. Default: max(10^3,10rho)  
**rho_factor:** factor by which rho changes each iteration. Default: 1.1  
**mu_factor:** factor by which mu changes each iteration. Default: 1.2  
**rho_max:** maximum value of rho. Default: 5 10^4  
**mu_max:** maximum value of mu. Default: 10^8  

Opts must be set. Algorithm options are:  
**QP_solver:** 'ipopt', 'linsolve', 'backslash', 'quadprog', 'pinv'. Pinv is best for non-convex problems and backslash is best for large-scale or convex problems  
**conv_Hess:** 1,0. Whether or not to convexify the Hessian matrices.  
**MIP_solver:** 'ipopt', 'bonmin', 'gurobi', 'cplex'  
**MINLPsolver_opts:** used to pass options to the MIP solver  
**NLPsolver_opts:** used to set ipopt options in re-solve phase  
**eig_flip:** 1,0. Whether or not to flip the eigenvalues of the Hessian matrices in the QP step  
**eig_non0:** 1,0. Whether or not to modify zero and negative eigenvalues (regularization) of the Hessian matrices in the QP step  

The syntax for using MI_ALADIN is shown in the EXAMPLE_MIP.

To cite use:  
@InProceedings{Murray_2018,  
  author    = {Murray, A. and Faulwasser, T. and Hagenmeyer, V.},  
  title     = {Mixed-Integer vs. Real-Valued Formulations of Distributed Battery Scheduling Problems},  
  booktitle = {10th Symposium on Control of Power and Energy Systems (CPES 2018)},  
  year      = {2018},  
  volume    = {51},  
  number    = {28},  
  pages     = {350--355},  
  }  
