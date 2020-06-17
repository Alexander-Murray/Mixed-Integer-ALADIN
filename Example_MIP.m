close all;
clear all;
clc;

import casadi.*

%% Problem parameters
N = 10; % number of partitions
NC = 5; % number of consensus constraints
nx = 10*ones(N,1); % number of real vars in each partition
nz = 10*ones(N,1); % number of discrete vars in each partition
nxz = nx+nz;

%% Problem setup
for i = 1:N
    x = SX.sym(strcat('x',num2str(i),'_'),nx(i),1); % local real vars
    z = SX.sym(strcat('z',num2str(i),'_'),nz(i),1); % local discrete vars
    
    % objective function
    obj{i} = rand*x.'*x + rand*z.'*z;
    obj_funs{i} = Function('f',{[x;z]},{obj{i}});
    
    % local constraints
    cons{i} = x + i*z - 5*rand;
    cons_funs{i} = Function('f',{[x;z]},{cons{i}});
    
    % consensus constraint matrices
    A_mats{i} = [rand(NC,nx(i)),zeros(NC,nz(i))];
    
    % initial guess
    x0{i} = zeros(nx(i)+nz(i),1);
    
    % box constraints
    lb{i} = -2*ones(nxz(i),1);
    ub{i} = 2*ones(nxz(i),1);
    
    % vector denoting which variables are integer valued
    discrete{i}=[zeros(1,nx(i)),ones(1,nz(i))];
end
b = zeros(NC,1); % RHS of consensus constraints
lam0 = zeros(NC,1); % initial guess of Lagrange multiplier of consensus constraints

%%
%Centralized solution
x_glob =  SX.sym('x',[sum(nxz),1]);

obj_cent_temp = 0;
cons_cent_temp= horzcat(A_mats{:})*x_glob-b;
for i = 1:N
    x_part = x_glob(sum(nxz(1:i-1))+1:sum(nxz(1:i)));

    obj_cent_temp = obj_cent_temp + obj_funs{i}(x_part);
    
    cons_cent_temp = [cons_cent_temp;cons_funs{i}(x_part)];
end
obj_cent = Function('f',{x_glob},{obj_cent_temp});

cons_cent = Function('g',{x_glob},{cons_cent_temp});

discrete_cent = horzcat(discrete{:});

lbg = [zeros(NC,1);-Inf(length(cons_cent_temp)-NC,1)];
ubg = [zeros(NC,1);zeros(length(cons_cent_temp)-NC,1)];

tic
% opts.verbose_init = 0;
% opts.bonmin.heuristic_RINS = 'yes';
% opts.bonmin.fp_pass_infeasible = 'yes';
% opts.bonmin.heuristic_feasibility_pump = 'yes';
% opts.bonmin.pump_for_minlp = 'yes';
% opts.bonmin.print_level = 0;
% opts.bonmin.bb_log_level = 0;
% opts.bonmin.fp_log_level = 0;
% opts.bonmin.lp_log_level = 0;
% opts.bonmin.milp_log_level = 0;
% opts.bonmin.nlp_log_level = 0;
% opts.bonmin.oa_log_level = 0;
% opts.bonmin.hessian_approximation = 'limited-memory';
% opts.bonmin.limited_memory_aug_solver = 'extended';
cent_opts.print_time = 0;
cent_opts.discrete = discrete_cent;

nlp =   struct('x',x_glob,'f',obj_cent(x_glob),'g',cons_cent(x_glob));
% cas =   nlpsol('solver','bonmin',nlp,opts);
cas =   qpsol('solver','gurobi',nlp,cent_opts);
sol =   cas('lbx', vertcat(lb{:}),...
            'ubx', vertcat(ub{:}),...
            'lbg', lbg,...
            'ubg', ubg,...
            'x0', vertcat(x0{:}));
Cent_time = toc;

cas.stats.return_status
Cent_sol = full(sol.x);

%% ALADIN solution
params.rho0 = 5*10^2;
params.max_iter = 100;
opts.QP_solver = 'backslash';
opts.MIP_solver = 'gurobi';
opts.conv_Hess = 1;
opts.eig_flip = 0;
opts.eig_non0 = 1;
params.rho_factor = 1.0;
params.mu_factor = 1.0;
opts.NLPsolver_opts.print_time = 0;
opts.MINLPsolver_opts.print_time = 0;
opts.NLPsolver_opts.ipopt.print_level = 0;
opts.NLPsolver_opts.ipopt.hessian_approximation = 'limited-memory';
opts.NLPsolver_opts.ipopt.limited_memory_aug_solver = 'extended';

[ ALADIN_xopt, logg ] = MI_ALADIN(obj_funs,cons_funs,A_mats,b,x0,lam0,lb,ub,discrete,params,opts);

ALADIN_time = sum(logg.par_time) + sum(logg.seq_time);

%% display solution
plot_length = size(logg.delY,2)+1;

% Plot Solution against optimal value and plot constraint violation
InOpt = zeros(1,plot_length);
ConViol = zeros(1,plot_length);
for i = 1:plot_length
   InOpt(i) = full(obj_cent(logg.X(:,i+1)));
   ConViol(i) = sum(abs(horzcat(A_mats{:})*logg.X(:,i+1)));
end

x = 1:plot_length;
close all;

figure(1);
y = InOpt(1:plot_length);
hold on;
grid on;
plot(x,y,plot_length,full(obj_cent(full(sol.x))),'-mo','MarkerFaceColor',[.49 0 .63]);
ylabel('Objective Function Value')
xlabel('Iterations')
xticks(0:ceil(plot_length/20):plot_length)
set(gcf,'position',[50 600 560 420])

figure(2);
z = ConViol(1:plot_length);
hold on;
grid on;
plot(x,z);
ylabel('$||\sum_iP_i(k)-G(k)||$','Interpreter','Latex')
xlabel('Iterations')
xticks(0:ceil(plot_length/20):plot_length)
set(gcf,'position',[600 600 560 420])


disp('Bonmin time')
disp(Cent_time)
disp('ALADIN time')
disp(ALADIN_time)
disp('ALADIN obj. / Bonmin obj.')
disp(full(obj_cent(ALADIN_xopt)/sol.f))
disp('|ALADIN x - Bonmin x|')
disp(sum(abs(ALADIN_xopt-full(sol.x))))
disp('Final ALADIN equality constraint violation')
disp(ConViol(plot_length))
disp('Min_epsilon')
disp(max(sum(abs(logg.X(plot_length)-logg.X(plot_length+1))),ConViol(plot_length)))

