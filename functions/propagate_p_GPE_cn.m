function [p_k] = propagate_p_GPE_cn(p_kp1, psi_k, psi_kp1, dt, V_k, V_kp1, m, g1d, mue, lap4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagates the adjunct state one step backwards
% 
% Args:
%     p_kp1 (array):                        adjunct state at k+1
%     psi_k (double):                       wavefunction at k
%     psi_kp1 (fun):                        wavefunction at k+1
%     V_k (double):                         potential at V(k)
%     V_kp1 (double):                       potential at V(k+1)
%     m (double):                           normalized mass
%     g1d (array):                          coupling constant g1D
%     mue (double):                         chemical potential
%     lap4 (array):                         N_x x N_x laplace matrix
% Returns:
%     p_k (array):                          adjunct state at k
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



n               = numel(V_k);

psi_k_q_real    = real(psi_k.^2);
psi_k_q_imag    = imag(psi_k.^2);
psi_kp1_q_real  = real(psi_kp1.^2);
psi_kp1_q_imag  = imag(psi_kp1.^2);

Ak_11   = speye(n) + dt/2*g1d*spdiag(psi_k_q_imag);
Ak_12   = dt/2*(-1/2/m * lap4 + spdiag(V_k - mue + 2*g1d*abs(psi_k).^2 - g1d*psi_k_q_real));
Ak_21   = dt/2*(1/2/m * lap4 - spdiag(V_k - mue + 2*g1d*abs(psi_k).^2 + g1d*psi_k_q_real));
Ak_22   = speye(n) - dt/2*g1d*spdiag(psi_k_q_imag);

Akp1_11 = speye(n) - dt/2*g1d*spdiag(psi_kp1_q_imag);
Akp1_12 = -dt/2*(-1/2/m * lap4 + spdiag(V_kp1 - mue + 2*g1d*abs(psi_kp1).^2 - g1d*psi_kp1_q_real));
Akp1_21 = -dt/2*(1/2/m * lap4 - spdiag(V_kp1 - mue + 2*g1d*abs(psi_kp1).^2 + g1d*psi_kp1_q_real));
Akp1_22 = speye(n) + dt/2*g1d*spdiag(psi_kp1_q_imag);

M_11    = Ak_11-Ak_12*spdiag((1./diag(Ak_22)))*Ak_21;
M_22    = Ak_22-Ak_21*spdiag((1./diag(Ak_11)))*Ak_12;
V_1     = Akp1_11 * real(p_kp1) + Akp1_12 * imag(p_kp1);
V_2     = Akp1_21 * real(p_kp1) + Akp1_22 * imag(p_kp1);

sol     = [M_11\V_1-M_11\(Ak_12*(spdiag(1./diag(Ak_22))*V_2));-M_22\(Ak_21*(spdiag(1./diag(Ak_11))*V_1))+M_22\V_2];
% sol = [M_11,zeros(n);zeros(n),M_22]\([speye(n),-Ak_12*spdiag(1./diag(Ak_22));-Ak_21*spdiag(diag(1./Ak_11)),speye(n)]*[V_1;V_2]);


p_k     = sol(1:n)+1i*sol(n+1:2*n);

end
