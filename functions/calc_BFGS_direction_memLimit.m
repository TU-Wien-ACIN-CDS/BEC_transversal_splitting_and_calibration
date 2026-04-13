function [p_k] = calc_BFGS_direction_memLimit(y_old,s_old,n_old,grad_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the search direction using the memory limited BFGS algorithm (see Wikipedia "Limited-memory BFGS")
% 
% Args:
%     ys_old (array):                       N_old_max x N_t array with changes in trajectory of previous steps (memory limited BFGS)
%     ss_old (array):                       N_old_max x N_t array with changes in gradient of previous steps (memory limited BFGS)
%     n_old (double):                       number of previous BFGS steps (can not be higher than N_old_max)
%     grad_k (array):                       N_t x 1 current gradient
%
% Returns:
%     p_k (array):                          N_t x 1 current search direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_old == 0
    p_k = -grad_k; % Standard steepest descent
    return;
end

q               = grad_k;
alpha           = zeros(1,n_old);
rho             = zeros(1,n_old);
for j = 1:n_old
    yj          = y_old(:,j);
    sj          = s_old(:,j);
    rhoj        = 1/(1e-15 + yj'*sj);
    alphaj      = rhoj*sj'*q;
    q           = q - alphaj*yj;
    rho(j)      = rhoj;
    alpha(j)    = alphaj;
end

s_newest = s_old(:, 1);
y_newest = y_old(:, 1);
y_norm_sq = y_newest' * y_newest;

if y_norm_sq > 1e-8
    gamma = (s_newest' * y_newest) / y_norm_sq;
else
    gamma = 1.0; % Fallback if gradient change is nearly zero
end
z               = gamma * q;
beta            = zeros(1,n_old);
for j = n_old:-1:1
    beta(j)     = rho(j)*y_old(:,j)'*z;
    z           = z + s_old(:,j)*(alpha(j) - beta(j));
end

if any(isnan(z))
    disp('search direction is NaN!!!!!!!!')
end

p_k = -z;

% if n_old == 0
    % p_k = -grad_k;
    % return;
% end
% 
% % Parameters
% tol_ys = 1e-12;
% 
% % Two-loop recursion (standard)
% q = grad_k;                    % start with negative gradient
% alpha = zeros(1, n_old);
% rho = zeros(1, n_old);
% 
% % Backward loop: most recent -> oldest (assume col 1 is most recent)
% for j = 1:n_old
%     yj = y_old(:, j);
%     sj = s_old(:, j);
%     ys = yj' * sj;
%     if ys <= tol_ys || all(sj==0) || all(yj==0)
%         rho(j) = 0;
%         alpha(j) = 0;
%         continue;               % skip pair with nonpositive curvature
%     end
%     rhoj = 1 / ys;
%     alphaj = rhoj * (sj' * q);
%     q = q - alphaj * yj;
%     rho(j) = rhoj;
%     alpha(j) = alphaj;
% end
% 
% % Choose H0 scalar gamma using the most recent valid pair (col 1)
% s_newest = s_old(:, 1);
% y_newest = y_old(:, 1);
% y_norm_sq = y_newest' * y_newest;
% if y_norm_sq > tol_ys
%     gamma = (s_newest' * y_newest) / y_norm_sq;
%     if gamma <= 0
%         gamma = 1.0;
%     end
% else
%     gamma = 1.0;
% end
% 
% z = gamma * q;
% 
% % Forward loop: oldest -> most recent
% for j = n_old:-1:1
%     if rho(j) == 0
%         continue;
%     end
%     beta = rho(j) * (y_old(:, j)' * z);
%     z = z + s_old(:, j) * (alpha(j) - beta);
% end
% 
% p_k = -z;
% end

% if n_old == 0
%     p_k = grad_k;
%     return
% end
% 
% alpha   = zeros(n_old,1);
% rho     = zeros(n_old,1);
% q       = grad_k;
% for i = 1:n_old
%     rho(i) = 1.0 / (y_old(:,i)'*s_old(:,i) + 1e-12);
%     alpha(i) = rho(i) * s_old(:,i)' * q;
%     q = q - alpha(i) * y_old(:,i);
% end
% gamma = (s_old(:,1)'* y_old(:,1) / (y_old(:,1)'* y_old(:,1) + 1e-12));
% r = gamma * q;
% for i = n_old:-1:1
%     beta = rho(i) * y_old(:,i)'* r;
%     r = r + s_old(:,i) * (alpha(i) - beta);
% end
% 
% p_k = -r;

% if n_old == 0
%     p_k = -grad_k; % Standard steepest descent
%     return
% end
% 
% q = grad_k;
% alpha = zeros(n_old, 1);
% rho = zeros(n_old, 1);
% 
% % LOOP 1: Newest to Oldest (1 -> n_old)
% for i = 1:n_old
%     si = s_old(:, i);
%     yi = y_old(:, i);
% 
%     rho(i) = 1.0 / (yi' * si + 1e-12);
%     alpha(i) = rho(i) * (si' * q);
%     q = q - alpha(i) * yi;
% end
% 
% % SCALING: Use the MOST RECENT step (Column 1)
% % This is the step most likely to represent the current local curvature.
% s_newest = s_old(:, 1);
% y_newest = y_old(:, 1);
% y_norm_sq = y_newest' * y_newest;
% 
% if y_norm_sq > 1e-12
%     gamma = (s_newest' * y_newest) / y_norm_sq;
% else
%     gamma = 1.0; % Fallback if gradient change is nearly zero
% end
% 
% r = gamma * q;
% 
% % LOOP 2: Oldest to Newest (n_old -> 1)
% for i = n_old:-1:1
%     si = s_old(:, i);
%     yi = y_old(:, i);
% 
%     beta = rho(i) * (yi' * r);
%     r = r + si * (alpha(i) - beta);
% end
% 
% p_k = -r;
% 
% p_k = - grad_k/1000;
% 
% end