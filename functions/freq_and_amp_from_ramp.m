function [slosh_y,psi_iter,A_t_full,t] = freq_and_amp_from_ramp(A_t,dt,V_d_A,x,par,u0,T_sloshing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates sloshing and breathing values for givenr amp
% 
% Args:
%     A_t (array):                          1 x N_t trajectory
%     dt (double):                          time discretization
%     V_d_A (fun):                          (A) calculates potential from A
%     x (array):                            spatial values
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
%     u0 (array):                           initial wave function
%     T_sloshing (double):                  Time used for sloshing
%
% Returns:
%     slosh_y (array):                      values of sloshing and breathing
%     psi_iter (array):                     N_x x N_t waveform values
%     A_t_full (array):                     1 x N_t parameter values including sloshing
%     t (array):                            1 x N_t time values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('T_sloshing','var')
    T_sloshing  = 2;
end

N_t_ramp                                    = numel(A_t);
t_ramp                                      = dt*(0:1:N_t_ramp-1);

t_sloshing                                  = t_ramp(end)+dt:dt:(t_ramp(end)+T_sloshing);
N_t_sloshing                                = numel(t_sloshing);
t                                           = [t_ramp,t_sloshing];
N_t_plot                                    = numel(t);
A_t_full                                    = [A_t,ones(1,N_t_sloshing)*A_t(end)];
psi_iter                                    = calc_prop_GPE(N_t_plot,x,u0,dt,V_d_A,A_t_full,par);

[relative_difference_iter,sigma_iter]       = calculate_difference_from_fit(x,psi_iter);
relative_diff_sloshing                      = relative_difference_iter(N_t_ramp+1:end);
sigma_sloshing                              = sigma_iter(1,N_t_ramp+1:end);

subsample_factor                            = 50;
fft_window_factor                           = 1e5;
relative_diff_sloshing_subsampled           = relative_diff_sloshing(1:subsample_factor:end)-mean(relative_diff_sloshing);
N_t_sloshing_subsampled                     = numel(relative_diff_sloshing_subsampled);
N_t_sloshing_subsampled_scaled              = 2^nextpow2(fft_window_factor*N_t_sloshing_subsampled);
Y_rel_diff                                  = fft(relative_diff_sloshing_subsampled,N_t_sloshing_subsampled_scaled)/N_t_sloshing_subsampled;
f_rel_diff                                  = 1/dt/subsample_factor/N_t_sloshing_subsampled_scaled*(0:N_t_sloshing_subsampled_scaled-1);

Y_half_diff                                 = 2*abs(Y_rel_diff(2:round(end/2)));
indices_local_max                           = islocalmax(Y_half_diff,2);
Y_half_diff_reduced                         = Y_half_diff;
Y_half_diff_reduced(~indices_local_max)     = 0;

[amp_maxk,maxk_indx_f]                      = maxk(Y_half_diff_reduced,2);
f_maxk                                      = f_rel_diff(maxk_indx_f+1);
if (amp_maxk(1)-amp_maxk(2))/amp_maxk(1) < 0.5
    [f_d,indx_small]                            = min(f_maxk);
    [f_d_2,~]                                   = max(f_maxk);
    amp_out                                     = amp_maxk(indx_small);
    phase_out                                   = mod(angle(Y_rel_diff(maxk_indx_f(indx_small)+1)),2*pi); % +pi/2 to shift from phase of cos to phase of sin
else
    f_d                                         = f_maxk(1);
    f_d_2                                       = f_maxk(2);
    amp_out                                     = amp_maxk(1);
    phase_out                                   = mod(angle(Y_rel_diff(maxk_indx_f(1)+1)),2*pi); % +pi/2 to shift from phase of cos to phase of sin
end
mean_d                                      = mean(relative_diff_sloshing);
var_d                                       = var(relative_diff_sloshing,0,'all');
mean_sigma                                  = mean(sigma_sloshing);


subsample_factor                            = 50;
fft_window_factor                           = 1e5;
breathing_subsampled                        = sigma_sloshing(1:subsample_factor:end)-mean(sigma_sloshing);
N_t_breathing_subsampled                    = numel(breathing_subsampled);
N_t_breathing_subsampled_scaled             = 2^nextpow2(fft_window_factor*N_t_breathing_subsampled);
Y_sigma                                     = fft(breathing_subsampled,N_t_breathing_subsampled_scaled)/N_t_breathing_subsampled;
f_sigma                                     = 1/dt/subsample_factor/N_t_breathing_subsampled_scaled*(0:N_t_breathing_subsampled_scaled-1);
% [amp,max_indx_f]                            = max(2*abs(Y_rel_diff(2:round(end/2))));

Y_sigma                                     = 2*abs(Y_sigma(2:round(end/2)));
indices_local_max                           = islocalmax(Y_sigma,2);
Y_half_sigma_reduced                        = Y_sigma;
Y_half_sigma_reduced(~indices_local_max)    = 0;

[amp_maxk,maxk_indx_f]                      = maxk(Y_half_sigma_reduced,2);
if (amp_maxk(1)-amp_maxk(2))/amp_maxk(1) < 0.5
    [max_indx_f,indx_is]                        = min(maxk_indx_f);
else
    max_indx_f                                  = maxk_indx_f(1);
    indx_is                                     = 1;
end

f_sigma                                     = f_sigma(max_indx_f+1);
amp_sigma                                   = amp_maxk(indx_is);

var_sigma                                   = var(sigma_sloshing,0,'all');
var_sigma_M                                 = var(1/sqrt(2)./sigma_sloshing,0,'all');

slosh_y                                     = [f_d;amp_out;mean_d;phase_out;var_d;mean_sigma;var_sigma;amp_sigma;f_sigma;var_sigma_M;f_d_2];

end