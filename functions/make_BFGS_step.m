function [lambda_t_new,ys_new,ss_new,n_new,gradients_new,search_directions,second_derivs,cost_init,a_opt,n_iter_opt,breaker,Hkp1] = make_BFGS_step(x,p_prop_function,psi_store,V_lambda,dV_lambda_dlambda,lambda_t,ys_old,ss_old,gradients_old,second_derivs_old,n_old,n_lambdas,cost_fun_barebones,par_bfgs,Hk,cost_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makes one BFGS line search step
% 
% Args:
%     x (array):                            N_x x 1 array with spatial values
%     p_prop_function:                      (psi,lambda) returns N_x x N_t array with adjunct state wavefunctions
%     psi_store (array):                    N_x x N_t array with state wavefunctions
%     V_lambda (fun):                       (lambda) gives the current potential as a Nx x 1 array
%     dV_lambda_dlambda (fun):              (lambda,n_lambda) gives the derivative of the potential for changes of the n_lambda-th lambda trajectory
%     lambda_t (array):                     N_lambda x N_t array with current trajectory
%     ys_old (array):                       N_old_max x N_t array with changes in trajectory of previous steps (memory limited BFGS)
%     ss_old (array):                       N_old_max x N_t array with changes in gradient of previous steps (memory limited BFGS)      
%     gradients_old (array):                N_old_max x N_t array with gradient of previous steps (memory limited BFGS)                   
%     second_deriv_old (array):             N_lambda x N_t array with old second derivatives of the trajectory
%     n_old (double):                       number of previous BFGS steps (can not be higher than N_old_max)
%     n_lambda (double):                    number of paramters for the potential
%     cost_fun_barebones (fun):             (lambda) calculates const-function value for given lambda value
%     par_bfgs (struct):                    BFGS parameters
%       - use_memLimit (flag):              (1 use memory limited BFGS, 0 use full BFGS (high RAM use for long trajectories))
%       - n_old_max (int):                  maximum number of old steps for memory limited BFGS
%       - gradient_second_deriv (fun):      (psi,p,dV_lambda,lambda,second_deriv_old) calculates the second derivative of the gradient
%       - lap_t (array):                    N_t x N_t array with laplace matrix in time
%       - use_numerical_derivative (flag):  decides if numerical derivative or dV_lambda_dlambda shall be used (1 use numerical, 0 use dV_lambda_dlambda)
%       - dt (double):                      Time discretization
%       - A_max (double):                   max value for trajectory
%       - A_min (double):                   min value of trajectory
%       - calc_second_deriv_t (fun):        (lambda) calculates second derivative in time of trajectories
%       - use_H1_cost (flag):               decides which gradient is used for(1 use H1 cost (see "Computational techniques for a quantum control problem with H1-cost"), 0 use L2 cost, 2 use combined H1 and L2 cost, 3 use smoothing of L2 cost)
%       - dt (double):                      time discretization
%       - C_H1_cost (double):               % scaling factor of L2 component of combined cost, only used if par_bfgs.use_H1_cost == 2
%     Hk (array):                           Hessian matrix approximation (only necessary for full BFGS)
%
% Returns:
%     lambda_t_new (array):                 N_lambda x N_t array new trajectories
%     ys_new (array):                       N_old_max x N_t array with changes in trajectory of previous and current steps (memory limited BFGS)
%     ss_new (array):                       N_old_max x N_t array with changes in gradient of previous and current steps (memory limited BFGS)      
%     gradients_new (array):                N_old_max x N_t array with gradient of previous and current steps (memory limited BFGS)          
%     search_directions (array):            N_t x N_lambda array with search directions of current step          
%     second_derivs (array):                N_t x N_lambda array with second derivs new trajectories
%     val_opt (double):                     current cost-function value 
%     a_opt (double):                       current step distance
%     n_iter_opt (int):                     number of optimizer iterations
%     breaker (flag):                       (1 breaks the outside loop, 0 does not break the outisde loop)
%     Hkp1 (array):                         new Hessian matrix approximation (only necessary for full BFGS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[N_t,~]             = size(gradients_old);
gradients_new       = gradients_old;

[p_store,psi_store_again]   = p_prop_function(psi_store,lambda_t); % calculates adjunct state for current trajectory

search_directions   = zeros(N_t,n_lambdas);
second_derivs_new   = second_derivs_old;

breaker             = 0;

% calculating search direction
for iter_lambdas = 1:n_lambdas
    y_old               = zeros(N_t,par_bfgs.n_old_max);
    y_old(:,:)          = ys_old(:,:,iter_lambdas);
    s_old               = zeros(N_t,par_bfgs.n_old_max);
    s_old(:,:)          = ss_old(:,:,iter_lambdas);
    second_deriv_old    = second_derivs_old(iter_lambdas,:);
    [grad_k,second_deriv]   = calc_Jlambda_grad(p_store,psi_store_again,V_lambda,dV_lambda_dlambda,par_bfgs.gradient_second_deriv,second_deriv_old,lambda_t,x,par_bfgs.dt,par_bfgs.use_numerical_derivative,n_lambdas,par_bfgs.use_H1_cost,par_bfgs.C_H1_cost);
    
    if ~all(gradients_old == 0)
        y_old(:,1)          = grad_k - gradients_old(:,iter_lambdas);
    end
    y_old([1,N_t],1)    = 0;
    s_old([1,N_t],1)    = 0;
    % y_old(:,1)          = y_old(:,1)/max(abs(y_old(:,1)),[],'all');
    ys_old(:,1,iter_lambdas)            = y_old(:,1);
    gradients_new(:,iter_lambdas)       = grad_k;
    if ~par_bfgs.use_memLimit
        [p_k,Hkp1]    = calc_BFGS_direction_full(y_old,s_old,Hk,grad_k);
    else
        p_k           = calc_BFGS_direction_memLimit(y_old,s_old,n_old,grad_k);
    end
    if (grad_k' * p_k) >= 0
        fprintf('Direction is NOT a descent direction! Resetting to steepest descent.\n');
        p_k     = -grad_k; 
        n_old   = 0; % Clear memory because it's corrupted
    end
    search_directions(:,iter_lambdas) = p_k;
    if any(isnan(search_directions(:,iter_lambdas)),'all')
        fprintf('search direction contains undefined values \n')
        search_directions(isnan(search_directions)) = 0;
    end
end


% defining maximum and minimum values for trajectory (should by set to +/- infinity if not used)
A_max           = par_bfgs.A_max;
A_min           = par_bfgs.A_min;

figure(11)
subplot(4,1,2)
plot(par_bfgs.dt*(1:N_t),search_directions)
grid on
hold on
title('search direction')
xlim([0,par_bfgs.dt*N_t])
xlabel('t in ms')
ylabel('delta A')
drawnow


% defining cost function for line search with step size a as an input 
cost_fun        = @(a) cost_fun_barebones(max(min(lambda_t + a.*search_directions',A_max),A_min));

a_init              = 1*ones(n_lambdas,1);
reduction_factor    = 0.8;
fprintf("a = 0.0000, cost = %.6f\n",cost_0)

for i = 0:50
    cost_init = cost_fun(a_init);
    fprintf("a = %.4f, cost = %.6f\n",a_init,cost_init)
    if cost_init < cost_0
        break
    end
    a_init = reduction_factor * a_init;
end
                    = a_test(indx_a_init);

% defines optimization options
% fminunc_opt     = optimoptions('fminunc','Display','iter-detailed');
fmincon_opt     = optimoptions('fmincon','Display','iter-detailed');
lb              = -eps;
ub              = 20;

fprintf("starting optimization with a = %.4f, and J = %.6f, (J0 = %.6f)\n",a_init,cost_init,cost_0)

% [a_opt,cost_init,outputflag,outputoptions] = fminunc(cost_fun,a_init,fminunc_opt);
[a_opt,cost_init,outputflag,outputoptions] = fmincon(cost_fun,a_init,[],[],[],[],lb,ub,[],fmincon_opt);

fprintf('optimal a: %f \n',a_opt)
n_iter_opt = outputoptions.iterations;

% calculates solution of current line search
lambda_t_new    = max(min(lambda_t + a_opt.*search_directions',A_max),A_min);

second_derivs   = second_derivs_old;
ss_new          = ss_old;
ys_new          = ys_old;

if abs(a_opt) > 1e-5
    
    n_new           = min(n_old + 1,par_bfgs.n_old_max);
    % upates arrays for BFGS updates
    for iter_lambdas = 1:n_lambdas
        s_new                    = (lambda_t_new(iter_lambdas,:)-lambda_t(iter_lambdas,:))';
        if (p_k'*s_new)<=0
            fprintf('curvature condition is not ensured!\n')
        end
        ss_new(:,:,iter_lambdas) = [s_new,ss_old(:,1:end-1,iter_lambdas)];
        ys_new(:,:,iter_lambdas) = [y_old(:,1),ys_old(:,1:end-1,iter_lambdas)];
        second_derivs(iter_lambdas,:) = par_bfgs.calc_second_deriv_t(lambda_t_new(iter_lambdas,:));
        % second_derivs(iter_lambdas,:) = second_derivs(iter_lambdas,:) - linspace(second_derivs(iter_lambdas,1),second_derivs(iter_lambdas,end),N_t);
    end

elseif abs(second_derivs_old(1,:) * grad_k) < 1e-5
    fprintf('clearing BFGS buffer because of small curvature!!\n')
    n_new = 0;
    for iter_lambdas = 1:n_lambdas
        s_new                    = (lambda_t_new(iter_lambdas,:)-lambda_t(iter_lambdas,:))';
        if (p_k'*s_new)<=0
            fprintf('curvature condition is not ensured!\n')
        end
        ss_new(:,:,iter_lambdas) = [s_new,ss_old(:,1:end-1,iter_lambdas)];
        ys_new(:,:,iter_lambdas) = [zeros(size(grad_k)),ys_old(:,1:end-1,iter_lambdas)];
        second_derivs(iter_lambdas,:) = smooth(second_derivs_old(iter_lambdas,:) + a_opt(iter_lambdas) * second_derivs_new(iter_lambdas,:),'rloess');
        % second_dervis(iter_lambdas,:) = second_derivs(iter_lambdas,:) - linspace(second_derivs(iter_lambdas,1),second_derivs(iter_lambdas,end),N_t);
    end
else
    fprintf('small step detected, not added to buffer!!\n')
    n_new           = n_old;
    second_derivs   = second_derivs_old;
    ss_new          = ss_old;
    ys_new          = ys_old;
end

if par_bfgs.use_memLimit
    Hkp1            = Hk;
else
    Hkp1            = Hk;
    n_new           = n_old;
end

end

