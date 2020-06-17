function [Sig,eps,maxiter,rho,mu,rho_factor,mu_factor,rho_max,mu_max] = MI_AL_param_check(params,x0,obj)
    if isfield(params,'Sig')
        Sig = params.Sig;
    else
        for i = 1:length(x0)
            Sig{i} = eye(length(x0{i}));
        end
    end
    if isfield(params,'eps')
        eps = params.eps;
    else eps = 10^-3;
    end
    if isfield(params,'maxiter')
        maxiter = params.maxiter;
    else maxiter = 300;
    end
    if isfield(params,'rho')
        rho = params.rho0;
    else rho = 0.1*obj(vertcat(x0{:}));
    end
    if isfield(params,'mu')
        mu = params.mu0;
    else mu = max(10*rho,10^3);
    end
    if isfield(params,'rho_factor')
        rho_factor = params.rho_factor;
    else rho_factor = 1.1;
    end
    if isfield(params,'mu_factor')
        mu_factor = params.mu_factor;
    else mu_factor = 1.2;
    end
    if isfield(params,'rho_max')
        rho_max = params.rho_max;
    else rho_max = 5*10^4;
    end
    if isfield(params,'mu_max')
        mu_max = params.mu_max;
    else mu_max = 1*10^8;
    end
end
