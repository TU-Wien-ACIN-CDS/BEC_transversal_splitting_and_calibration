function [p_k,Hkp1] = calc_BFGS_direction_full(y_old,s_old,Hk,grad_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the search direction using the full BFGS algorithm
% 
% Args:
%     y_old (array):                        1 x N_t array with changes in trajectory of the previous step
%     s_old (array):                        1 x N_t array with changes in gradient of the previous step
%     H_k (array):                          N_t x N_t previous Hessian matrix aproximation
%     grad_k (array):                       N_t x 1 current gradient
%
% Returns:
%     p_k (array):                          N_t x 1 current search direction
%     Hkp1 (array):                         N_t x N_t new Hessian matrix aproximation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sk      = s_old(:,1);
yk      = y_old(:,1);

rhok    = 1/(1e-15 + yk'*sk);
n       = numel(grad_k);
Hkp1    = (eye(n)-rhok*sk*yk')*Hk*(eye(n)-rhok*yk*sk')+rhok*(sk*sk');

z       = Hkp1*grad_k;

p_k     = -z/max(abs(z),[],'all');

subplot(5,1,5)
plot(p_k)
grid on
hold on
title('search direction')

end