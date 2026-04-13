function J = cost_val_multi(cost_fun,N_eta,x,u0,dt,V_d,eta_t,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the cost-function value
% 
% Args:
%     cost_fun (fun):                       (psi,lambda,t) calculates cost-function value
%     N_eta (int):                          Number of time points
%     x (array):                            spatial values
%     u0 (array):                           initial wave function
%     dt (double):                          time discretization
%     V_d (fun):                            (lamabda) calculates potential
%     eta_t (array):                        1 x N_t trajectory
%     par (struct):                         model parameters
%       - L (double):                       size of spatial dimension
%       - N_x (int):                        Number of spatial points
%       - dx (double):                      spatial discretization
%       - x (array):                        1 x N_x spatial values
%       - m (double):                       normalized mass
%       - ilap (array):                     1 x N_x approximated fourier transform of laplace operator
%       - g (double):                       coupling constant g1D
%       - gamma (double):                   gamma for cost function
%       - eta (double):                     second costfunction parameter
%       - dt (double):                      time discretization
%       - lap4 (array):                     N_x x N_x laplace operator in space
%       - ilap4 (array):                    1 x N_x fourth order approximated fourier transform of laplace operator
%     nonpoly (flag):                       (1 use nonpolynomial SE, 0 use GPE)
%     solver_vers (flag):                   (1 use GPE, 0 use crank-nichols newton, 2 use crank-nichols fsolve)
%
% Returns:
%     J (double):                           cost-function value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[psi_store_gradient] = calc_prop_GPE(N_eta,x,u0,dt,V_d,eta_t,par);

J = cost_fun(psi_store_gradient(:,end),eta_t,(0:N_eta-1)*dt);

if J < 0
    toc
end

end