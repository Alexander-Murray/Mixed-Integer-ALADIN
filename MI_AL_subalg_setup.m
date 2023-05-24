function [minlp,nlp] = MI_AL_subalg_setup(solver_struct,opts,discrete)
    import casadi.*
    opts.MINLPsolver_opts.discrete = discrete;
        if strcmp(opts.MIP_solver,'bonmin')
            minlp = nlpsol('solver','bonmin',solver_struct,opts.MINLPsolver_opts);
        elseif strcmp(opts.MIP_solver,'gurobi')
            minlp = qpsol('solver','gurobi',solver_struct,opts.MINLPsolver_opts);
        elseif strcmp(opts.MIP_solver,'cplex')
            minlp = qpsol('solver','cplex',solver_struct,opts.MINLPsolver_opts);
        elseif strcmp(opts.MIP_solver,'ipopt')
            if sum(discrete)~=0
               error('IPOPT can only be applied to real-valued problems. Non-zero discrete vector detected.') 
            end
            minlp = nlpsol('solver','ipopt',solver_struct,opts.MINLPsolver_opts);
        end
     if isfield(opts,'NLPsolver_opts')
        nlp = nlpsol('solver','ipopt',solver_struct,opts.NLPsolver_opts); %needed for obtaining langrange multipliers from MIP subprobs
     else 
        nlp = nlpsol('solver','ipopt',solver_struct); %needed for obtaining langrange multipliers from MIP subprobs
     end
end