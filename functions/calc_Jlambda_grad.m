function [Jlambda_grad,second_deriv] = calc_Jlambda_grad(p_store,psi_store,V_lambda,dV_lambda_dlambda,gradient_second_deriv,second_deriv_old,lambda_t,x,dt,use_numerical_derivative,n_lambda,use_H1_cost,C_H1_cost)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the gradient of the costfunction in regards to the trajectory lambda
% 
% Args:
%     p_store (array):                  N_x x N_t array with adjunct state wavefunctions
%     psi_store (array):                N_x x N_t array with state wavefunctions
%     V_lambda (fun):                   (lambda) gives the current potential as a Nx x 1 array
%     dV_lambda_dlambda (fun):          (lambda,n_lambda) gives the derivative of the potential for changes of the n_lambda-th lambda trajectory
%     gradient_second_deriv (fun):      (psi,p,dV_lambda,lambda,second_deriv_old) calculates the second derivative of the gradient
%     second_deriv_old (array):         N_lambda x N_t array with old second derivatives of the trajectory
%     lambda_t (array):                 N_lambda x N_t array with current trajectory
%     x (array):                        N_x x 1 array with spatial values
%     dt (double):                      Time discretization
%     use_numerical_derivative (flag):  decides if numerical derivative or dV_lambda_dlambda shall be used (1 use numerical, 0 use dV_lambda_dlambda)
%     n_lambda (double):                number of paramters for the potential
%     use_H1_cost (flag):               decides which function space is used for gradient calculation (0: L2 cost with forced boundary values, 1: H1 cost, 2: combined cost with C_H1_cost*gradient - second derivative of gradient = RHS, 3: L2 cost with smoothing to conform to boundary values
%     C_H1_cost (double):               Scaling factor for combined scaling factor of L2 component of combined cost, only used if par_bfgs.use_H1_cost == 2
%
% Returns:
%     Jlambda_grad (array):             N_lambda x N_t array with calculated gradient values
%     second_deriv (array):             N_lambda x N_t array with second derivative of gradient values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N_x         = numel(x);
N_t         = size(lambda_t,2);

% calculates derivatives of V with regards to lambda
dV_lambda   = zeros(N_x,N_t);
if use_numerical_derivative
    eps = 1e-6;
    for j = 1:N_t
        lambda_plus = lambda_t(:,j);
        lambda_plus(n_lambda) = lambda_plus(n_lambda) + eps;
        lambda_minus = lambda_t(:,j);
        lambda_minus(n_lambda) = lambda_minus(n_lambda) - eps;
        dV_lambda(:,j)  = (V_lambda(lambda_plus) - V_lambda(lambda_minus)) / (2*eps);
    end
else
    for j = 1:N_t
        dV_lambda(:,j)  = dV_lambda_dlambda(lambda_t(j),n_lambda);
    end
end

% calculates second derivative of gradient of J to lambda
[second_deriv,Jl_pp_wo]   = gradient_second_deriv(psi_store,p_store,dV_lambda,lambda_t(n_lambda,:),second_deriv_old);

switch use_H1_cost
    case 0 % L2 cost
        Jlambda_grad                = second_deriv * dt^2;
    case 1 % H1 cost
        Lap_t                       = (spdiags(1,-1,N_t,N_t)+spdiags(-2,0,N_t,N_t)+spdiags(1,1,N_t,N_t));
        % Lap_t(1,1:3)                = [1,-2,1];
        % Lap_t(end,end-2:end)        = [1,-2,1];
        Jlambda_grad                = - Lap_t \ second_deriv * dt^2  * dt;
    case 2 % H1 cost extended with L2 cost (still type of H1 cost actually)
        Lap_t                       = (spdiags(1,-1,N_t,N_t)+spdiags(-2,0,N_t,N_t)+spdiags(1,1,N_t,N_t));
        Jlambda_grad                = (C_H1_cost*speye(N_t) - Lap_t) \ second_deriv * dt^2 * dt;
    case 3 % L2 cost with smoothing window to confirm with boundary values
        window_sigma                = 50;
        time_indices                = 1:N_t;
        smoothing_window            = erf(time_indices/(sqrt(2)*window_sigma)) - erf((time_indices - N_t)/(sqrt(2)*window_sigma)) - 1;
        Jlambda_grad               = second_deriv .* smoothing_window' * dt^2;
end

Jlambda_grad                = Jlambda_grad - linspace(Jlambda_grad(1),Jlambda_grad(end),N_t)'; % subtracting linear component to conform with boundary values

% figure(13)
% subplot(5,1,1)
% plot(Jl_pp_wo)
% grid on
% hold on
% title('second derivative wo App')
% subplot(5,1,2)
% plot(second_deriv)
% grid on
% hold on
% title('second derivative')
% subplot(5,1,3)
% plot(Jlambda_grad_1)
% grid on
% hold on
% title('gradient')
% subplot(5,1,4)
% plot(Jlambda_grad)
% grid on
% hold on
% title('gradient with BV')

end