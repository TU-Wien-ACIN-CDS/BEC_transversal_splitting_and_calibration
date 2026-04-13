function [p_store, psi_p] = calc_prop_GPE_p_v2(psi_p, lambda_p, x, lap4, V, m, g1d, mue, t_p, end_val_fun, x_selector)

N_t                     = numel(t_p);
p_store                 = zeros(numel(x),N_t);

p_end                   = end_val_fun(psi_p(:,end),lambda_p(:,end));
p_store(:,end)          = p_end(:);

% mue = g1d * max(abs(p_end)).^2;
% u_norm  = @(u)u/sqrt(trapz(x,conj(u).*u));

for j = N_t:-1:2
    p                       = propagate_p_GPE_cn(p_store(:,j),psi_p(:,j-1),psi_p(:,j),t_p(j)-t_p(j-1),V(lambda_p(j-1)),V(lambda_p(j)),m,g1d,mue,lap4);
    % p                       = u_norm(p);
    p_store(:,j-1)          = p(:);
end





end