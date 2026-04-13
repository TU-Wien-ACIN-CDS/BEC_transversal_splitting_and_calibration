function [J_lambda_pp,J_lambda_pp_without_lambda_pp] = calculate_J_lambda_pp(psi_store,p_store,dV_lambda,lambda,second_deriv_old,x,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate second derivative of cost function gradient
% 
% Args:
%     psi_store (array):                    N_t x N_t wave functions of psi
%     p_store (array):                      N_t x N_t wave functions of adjunct state p
%     dV_lambda (array):                    N_t x N_t derivative of Potential to trajcetroy values
%     lambda (array):                       1 x N_t trajectory
%     second_deriv_old (array):             1 x N_t second derivative of trajectory
%     x (array):                            1 x N_x spatial values
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
%
% Returns:
%     J_lambda_pp (array):                      N_t x 1 second derivative of J gradient
%     J_lambda_pp_without_lambda_pp (array):    N_t x 1 second derivative of J gradient without regularization from current trajectory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrand_term = real(conj(p_store).*dV_lambda.*psi_store);
J_lambda_pp_without_lambda_pp   = trapz(x,integrand_term,1)';
J_lambda_pp                     = -2*par.gamma*second_deriv_old' - J_lambda_pp_without_lambda_pp;


figure
subplot(2,3,1)
imagesc(1:numel(lambda),x,abs(p_store).^2)
subplot(2,3,2)
imagesc(1:numel(lambda),x,abs(psi_store).^2)
subplot(3,3,3)
imagesc(1:numel(lambda),x,dV_lambda)
subplot(2,3,4)
imagesc(1:numel(lambda),x,angle(p_store))
subplot(2,3,5)
imagesc(1:numel(lambda),x,angle(psi_store))
subplot(3,3,6)
imagesc(1:numel(lambda),x,real(conj(p_store).*psi_store))
subplot(3,3,9)
imagesc(1:numel(lambda),x,integrand_term)

end