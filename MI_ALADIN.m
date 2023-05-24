function [ xopt, log , tot_time ] = MI_ALADIN( obj_funs,cons_funs,A_mats,b,x0,lam0,lbx,ubx,discrete,params,opts)

good_statuses = {'Solve_Succeeded','integer optimal, tolerance','OPTIMAL','SUCCESS','integer optimal solution', 'Solved_To_Acceptable_Level','Maximum_Iterations_Exceeded'};                                

% get dimensions
NsubSys = length(obj_funs); %number of partitions
Ncons   = size(A_mats{1},1); %number of equality constraints

% build global A matrix
A       = horzcat(A_mats{:});

% build centralized objective function
xx        = SX.sym('x', [length(vertcat(x0{:})), 1]);
cent_obj  = 0;
var_count = 0;
nx = zeros(NsubSys,1);
for fun = 1:NsubSys
    nx(fun)   = length(x0{fun}); %number of variables in partition i
    cent_obj  = cent_obj + obj_funs{fun}(xx(var_count+1:var_count+nx(fun)));
    var_count = var_count+nx(fun);
end  
cent_obj_fun = Function('f',{xx},{cent_obj});

% set up params
[Sig,eps,maxiter,rho,mu,rho_factor,mu_factor,rho_max,mu_max] = MI_AL_param_check(params,x0,cent_obj_fun);
QPsolver = opts.QP_solver;

%% build local subproblems and CasADi functions
import casadi.*
rhoCas = SX.sym('rho',1,1);
lamCas = SX.sym('lam',length(lam0),1);

Hess_Loc_Fun = cell(1,NsubSys);
tic
for i=1:NsubSys
    xCas     = SX.sym('y',nx(i),1); %local decision variables
    diff_x   = find(ones(1,nx(i))-discrete{i}); %identify real variables for differentiation
    yCas     = SX.sym('x',nx(i),1); %parameterize QP solution
    x_opt{i} = x0{i}; %initial point given to local solvers
    
    % parameter vector of CasADi
    pCas = [ rhoCas;
             lamCas;
             yCas];
                
    % objective function for local NLP's
    obj_fun_Loc_Cas = obj_funs{i}(xCas) + lamCas'*A_mats{i}*xCas ...
                + rhoCas/2*(xCas - yCas)'*Sig{i}*(xCas - yCas);

    
    %Gradient of local objective
    gCas     = gradient(obj_funs{i}(xCas),xCas(diff_x));
    g_fun{i} = Function(['g' num2str(i)],{xCas},{gCas});
    
    
    % local inequality constraints
    cons_Cas  = cons_funs{i}(xCas);

    % Jacobian of local constraints
    Jac_Cas      = jacobian(cons_Cas,xCas(diff_x));   
    Jac_Fun{i}   = Function(['Jac' num2str(i)],{xCas},{Jac_Cas});
    Nloc_cons{i} = size(Jac_Cas,1); %number of local inequality constraints

    %Hessian  of local objective 
    kappa_Cas       = SX.sym('kapp',Nloc_cons{i},1);
    rho_Cas         = SX.sym('rho',1,1);
    Hess_Loc_Cass   = hessian(obj_funs{i}(xCas)+kappa_Cas'*cons_Cas,xCas(diff_x)) + opts.conv_Hess*rho_Cas*Sig{i}(diff_x,diff_x);     %The third term of the Hessian is to ensure positive-definiteness. Otherwise, the QP will be non-convex
    Hess_Loc_Fun{i} = Function(['H' num2str(i)],{xCas,kappa_Cas,rho_Cas},{Hess_Loc_Cass}); 

    % set up local solvers
    solver_struct = struct('x',xCas,'f',obj_fun_Loc_Cas,'g',cons_Cas,'p',pCas);
    [minlp{i},nlp{i}] = MI_AL_subalg_setup(solver_struct,opts,discrete{i});    
end
toc

%% Initialization
i        = 1; % iteration number
yy       = x0; % primal solution
lam      = lam0; % dual solution
delx     = inf; % QP solution
tot_time = 0; % total runtime
infeas   = 0; % problem feasibility token

%track algorithm progress
ObjVal = full(cent_obj_fun(vertcat(x0{:})));
ConViol = full(A*vertcat(yy{:}) - b); 
MI_AL_plot(0,ObjVal,ConViol,maxiter)
 
%setup a log for some key vectors
log         = struct();
log.X       = vertcat(x0{:});
log.delY    = [];   
log.Kappa   = [];
log.lambda  = lam;
log.ObjVal  = full(cent_obj_fun(log.X(:,1)));
log.ConViol = ConViol;
log.status  = 'running'; %this status should not appear
log.par_time = [];
log.seq_time = [];

%%

while i <= maxiter  
    
    sub_time = zeros(NsubSys,1);

    for j=1:NsubSys
        tic
        % set up parameter vector for local NLP's
        pNum = [ rho;
                 lam;
                 yy{j}]; 
             
        % solve local MINLP's   
        [x_opt{j},kappa_opt{j},status] = MI_AL_MINLP(minlp{j},nlp{j},x_opt{j},pNum,lbx{j},ubx{j},discrete{j},Nloc_cons{j},opts.MIP_solver);
        if ~any(strcmp(good_statuses,status))
            if i>1
               warning('rho_max may be set too low or rho_update too high.') 
            end
            infeas = 1;
            keyboard
            warning(  '%s, Iteration: %i, Partition: %i',minlp{j}.stats.return_status,i,j)
            break
        end
      
        % active set detection
        subprob_cons = full(cons_funs{j}(x_opt{j}));
        inact = subprob_cons<-1e-5;
        kappa_opt{j}(inact)= 0;

        % Jacobians of inequality constraints (lower left component of AQP)
        Ci = full(Jac_Fun{j}(x_opt{j}));
        for d = 1:length(discrete{j})
            if discrete{j}(d)==1
                Ci(:,d) = 0;
            end
        end
        C{j} = Ci(~inact,:); % eliminate inactive entries

        % evaluate gradients and hessians for the QP
        Hess_Loc_FunEval{j} = Hess_Loc_Fun{j}(x_opt{j},kappa_opt{j},rho);
        g_fun_eval{j}       = g_fun{j}(x_opt{j});
        
        if or(opts.eig_flip,opts.eig_non0)
            % eigenvalue decomposition of the hessian
            [V,D]       = eig(full(Hess_Loc_FunEval{j}));
            e           = diag(D);
            if opts.eig_non0
                % modify zero and negative eigenvalues (regularization)
                reg         = 1e-4;
                e(e<reg)    = reg; % 1e-4 
            else
                % flip the eigenvalues
                e           = abs(e); 
            end
            Hess_Loc_FunEval{j}  = V*diag(e)*transpose(V);
        end
        sub_time(j) = toc;
    end       
    tot_time = tot_time + max(sub_time);
    log.par_time = [log.par_time, sub_time];
    if infeas == 1
       log.status = 'INFEASIBLE';
       break
    end
    
    %Termination Check (only consensus, not objective value change)
    if max(abs(horzcat(A_mats{:})*vertcat(x_opt{:})-b))<eps
        xopt = vertcat(x_opt{:});
        log.status  = 'OPTIMAL';
        log.X          = [log.X vertcat(x_opt{:})];
        log.Kappa      = [log.Kappa vertcat(kappa_opt{:})];
        log.ObjVal     = [log.ObjVal full(cent_obj_fun(vertcat(x_opt{:})))];
        log.ConViol    = [log.ConViol A*vertcat(x_opt{:}) - b];
       break 
    end
    
    % build the QP
    rhsQP  = b-A*vertcat(x_opt{:});    
    
    C_QP   = rref(blkdiag(C{:})); 
    
    %linsolve and pinv don't seem to allow for sparse matrices
    if or(strcmp(QPsolver,'linsolve'),strcmp(QPsolver,'pinv'))
        HQP     = full(blkdiag(Hess_Loc_FunEval{:},mu*eye(Ncons)));  
        gQP     = full(vertcat(g_fun_eval{:},lam));
    else
        HQP     = sparse(blkdiag(Hess_Loc_FunEval{:},mu*eye(Ncons)));
        gQP     = sparse(vertcat(g_fun_eval{:},lam));
    end    
        
    AQP     = [A  -eye(Ncons); C_QP  zeros(size(C_QP,1),Ncons)];
    bQP     = [rhsQP;zeros(size(C_QP,1),1)];
    
    %remove ints from AQP
    ints = find(horzcat(discrete{:}));
    AQP(:,ints)=[];
    
    %remove zero rows from AQP
    remove = find(sum(abs(AQP'))==0);
    AQP(remove,:)=[];
    bQP(remove)=[];
    
    % solve global QP
    tic
    [QP_sol, lamges] = solveQP(HQP,gQP,AQP,bQP,QPsolver);
    qp_time = toc;
    tot_time = tot_time + qp_time;
    log.seq_time = [log.seq_time, qp_time];
    delx = QP_sol(1:(end-Ncons)); 
    
alpha_1 = 1; alpha_2 = 1; alpha_3 = 1; % default values. Change these with line search

%% update vars and logging
    % lambda only for cons. constr.
    lam = lam + alpha_3*(lamges(1:Ncons) - lam);

    % update the real-valued local variables
    ctr = 1;
    for j=1:NsubSys
        delta_x{j} = zeros(nx(j),1);
        for var = 1:nx(j)
            if discrete{j}(var)==0
                delta_x{j}(var) = delx(ctr);
                ctr = ctr + 1;
            end
        end
        yy{j} = yy{j} + alpha_1*(x_opt{j}-yy{j}) + alpha_2*delta_x{j};  
    end 

    % logging
    log.X          = [log.X vertcat(x_opt{:})];
    log.delY       = [log.delY vertcat(delta_x{:})];
    log.Kappa      = [log.Kappa vertcat(kappa_opt{:})];
    log.lambda     = [log.lambda lam];
    log.ObjVal     = [log.ObjVal full(cent_obj_fun(vertcat(x_opt{:})))];
    log.ConViol    = [log.ConViol A*vertcat(x_opt{:}) - b];
    
    % plotting
    ObjVal = full(cent_obj_fun(log.X(:,i+1))); %cent_obj_fun(log.X(:,i+1)+log.delY(:,i));
    ConViol = full(max(abs(A*vertcat(x_opt{:}) - b)));
    MI_AL_plot(i,ObjVal,ConViol,maxiter)
    
    %next iteration
    i = i+1;
    
    %updating rho and mu slightly at each iteration can help improve
    %convergence, although this is a heuristic tool
    if rho < rho_max
       rho = rho*rho_factor;
       sprintf('rho: %d',rho)
    end
    if mu < mu_max
       mu = mu*mu_factor;
       sprintf('mu: %d',mu)
    end
    

end

if i>= maxiter
   log.status  = 'timeout'; 
end

close all;
xopt = vertcat(x_opt{:});


end