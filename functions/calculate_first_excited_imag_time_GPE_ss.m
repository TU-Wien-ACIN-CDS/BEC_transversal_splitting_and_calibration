function [u,i] = calculate_first_excited_imag_time_GPE_ss(V,u0,par,u1_est,x,dt_gs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates first excited state of potential V using imaginary step method
% and orthoginality condition
% 
% Args:
%     V (array):                            N_x x 1 Potential values
%     u0 (array):                           N_x x 1 ground state
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
%     u1_est (array):                       N_x x 1 initial guess of excited state
%     x (array):                            N_x x 1 spatial values
%     dt_gs (double):                       time discretization for imaginary step
%
% Returns:
%     u (array):                            N_x x 1 calculated first excited state
%     i (int):                              number of iterations taken
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g1d     = par.g;
m       = par.m;

N_max   = 1e6;
u       = -sign(x).*abs(u0);
u_orth  = @(u) u - trapz(x,conj(u0).*u)*u0;
u_norm  = @(u)u/sqrt(trapz(x,conj(u).*u));

ilap    = par.ilap;

if any(isnan(V))
    fprintf('NaN in V for GS')
end

if ~exist("dt_gs")
    dt_gs = 1e-5;
end

for i = 1:N_max
    u_next = propagate_GPE_ss(u, -1i*dt_gs,V,m,g1d,0,ilap);
    u_next = u_orth(u_next);
    u_next = u_norm(u_next);
    if trapz(x,abs(u_next - u))<1e-11
        fprintf('Excited state converged in %i steps\n',i)
        break;
    end
    if any(isnan(u_next))
        tic
    end
    u = u_next;
end

[~,indx_max] = max(abs(u));
u = u * exp(-1i*angle(u(indx_max)));

if i == N_max
    fprintf('N_max reached in excited calc\n')
end


end