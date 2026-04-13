function [psi_store] = calc_prop_GPE(N_eta,x,u0,dt,V_d,eta_t,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates psi using 
% 
% Args:
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
%     psi_store_gradient (array):           N_x x N_t wavefunctions of psi calculated using GPE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_x             = numel(x);
m               = par.m;
ilap_scaled     = par.ilap;

psi_store       = zeros(N_x,N_eta);
psi_store(:,1)  = u0;
u               = u0;
g               = par.g;
mue_0           = g*max(abs(u0).^2);

for i_t = 1:N_eta-1
    u = propagate_GPE_ss(u,dt,V_d(eta_t(:,i_t)),m,g,mue_0,ilap_scaled);
    u = u./sqrt(trapz(x,abs(u).^2));
    psi_store(:,i_t+1) = u;
end


end