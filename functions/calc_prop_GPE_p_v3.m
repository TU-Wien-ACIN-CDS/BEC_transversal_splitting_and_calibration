function [p_store,psi_p_internal] = calc_prop_GPE_p_v3(psi_p, lambda_p, x, lap_scaled, ilap_scaled, V, g1d, mue, t_p, end_val_fun)

N_t                     = numel(t_p);
p_store                 = zeros(numel(x),N_t);
psi_p_internal          = zeros(numel(x),N_t);

u_norm  = @(u)u/sqrt(trapz(x,conj(u).*u));
    
p_end                   = end_val_fun(psi_p(:,end),lambda_p(:,end));
% p_end(:)                = u_norm(p_end(:));
p_store(:,end)          = p_end(:);
psi_p_internal(:,end)   = psi_p(:,end);

conditions                  = zeros(1,N_t);

% mue = g1d * max(abs(p_end)).^2;

% mue = mue * 0;

n_sub           = 5;

p = p_store(:,end);
psi = psi_p_internal(:,end);

dt = (t_p(end-1)-t_p(end));

for j = N_t:-1:2
    dt_sub = dt/n_sub;
    for i = 1:n_sub
        lambda_sub              = 1/2 * (lambda_p(j-1) + lambda_p(j));
        [p,psi,conditions(j)]   = propagate_p_GPE_both(p,psi,dt_sub,V(lambda_sub)-mue,V(lambda_sub)-mue,g1d,lap_scaled,ilap_scaled);
        % [p,conditions(j)]       = propagate_p_GPE_midpoint(p_store(:,j),(psi_p(:,j-1)+psi_p(:,j))/2,t_p(j)-t_p(j-1),(V(lambda_p(j-1))-mue+V(lambda_p(j))-mue)/2,g1d,lap_scaled);
        % p                       = u_norm(p);
        % psi                     = u_norm(psi);
    end
    p_store(:,j-1)          = p(:);
    psi_p_internal(:,j-1)   = psi(:);
end
% 
% figure
% plot(conditions)
% grid on
% title('condition')



end