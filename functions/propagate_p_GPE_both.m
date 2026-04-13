function [p_k,psi_k,condition] = propagate_p_GPE_both(p_kp1, psi_kp1, dt, Vmue_k, Vmue_kp1, g1d, lap_scaled, ilap_scaled)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagates the adjunct state and the state one step backwards
% 
% Args:
%     p_kp1 (array):                        adjunct state at k+1
%     psi_kp1 (fun):                        wavefunction at k+1
%     V_k (double):                         potential at V(k)
%     V_kp1 (double):                       potential at V(k+1)
%     m (double):                           normalized mass
%     g1d (array):                          coupling constant g1D
%     mue (double):                         chemical potential
%     lap4 (array):                         N_x x N_x laplace matrix
% Returns:
%     p_k (array):                          adjunct state at k
%     psi_k (array):                        state at k
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%% copied from paper
%  decompose adjoint variable into real and imaginary part

ut = exp(-1i * ilap_scaled * dt);                           % (6)
uv = exp(-1i * (Vmue_k + g1d*abs(psi_kp1).^2) * dt/2);    % (8) + (10)

p_i = propagate(psi_kp1,p_kp1,Vmue_k,g1d,dt);             % (21a)
psi_i = uv .* psi_kp1;                                      % (22a)

p_ii = ifft(ut .* fft(p_i));                                % (21b)
psi_ii = ifft(ut .* fft(psi_i));                            % (22b)

uv = exp(-1i*(Vmue_k + g1d*abs(psi_ii).^2)*dt/2);         % (25)

psi_k = uv .* psi_ii;                                       % (27)
p_k = propagate(psi_k,p_ii,Vmue_k,g1d,dt);                % (28)

%%

if ~exist("condition",'var')
    condition = 0;
end



function  pt = propagate( psi, p, v, kappa, dt )
%  Borzi and Hohenester, SIAM J. Sci. Comp. 30, 441 (2008), Eqs. (22,23)

%  decompose adjoint variable into real and imaginary part
pr = real( p ); 
pi = imag( p );
%  time propagation vector
u = 0.5 * dt * [ 1i * kappa * real( psi   .^ 2 ),     ...
                  2 * kappa * abs(  psi ) .^ 2 + v,   ...
               - 1i * kappa * imag( psi   .^ 2 ) ];
%  length of vector
a = sqrt( sum( u .^ 2, 2 ) );
%  unit vector
e = u ./ repmat( a, 1, 3 );

%  update p
prp = ( cos( a ) + 1i * sin( a ) .* e( :, 3 ) ) .* pr +  ...
     sin( a ) .* ( 1i * e( :, 1 ) + e( :, 2 ) ) .* pi;
pip = ( cos( a ) - 1i * sin( a ) .* e( :, 3 ) ) .* pi +  ...
     sin( a ) .* ( 1i * e( :, 1 ) - e( :, 2 ) ) .* pr;

pt = prp + 1i * pip;


end

end