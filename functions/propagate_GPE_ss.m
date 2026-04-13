function [u] = propagate_GPE_ss(u, dt,V,m,g1d,mue,ilap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagates the GPE one step
% 
% Args:
%     u (array):                            previous wavefunction
%     dt (double):                          time discretization
%     V (fun):                              (A) calculates potential from A
%     m (double):                           normalized mass
%     g1D (double):                         coupling constant g1D
%     mue (double):                         chemical potential
%     ilap (array):                         1 x N_x approximated fourier transform of laplace operator
%
%
% Returns:
%     u (array):                            new wavefunction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = exp(-1i * (V + g1d*abs(u).^2-mue) * dt / 2) .* u;

u = fftn(u);

kin_en  = ( - 0.5 * ilap / m ); 
u       = exp( - 1i * kin_en * dt  ).* u;
    
u = ifftn(u);

u = exp(-1i * (V + g1d*abs(u).^2-mue) * dt / 2) .* u;

end