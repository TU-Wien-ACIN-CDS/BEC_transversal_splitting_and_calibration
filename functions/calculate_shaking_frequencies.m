function [sim_freq,sim_amp,sim_freq_x0,sim_freq_x0_both,sort_indx,last_out_struct,f_fft_x0,Y_fft_x0,sim_freq_x0_both_maxk] = calculate_shaking_frequencies(As,V_d_A,x,par,N_cases,u0_list,dt,indx_offset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates frequencies from an offset/shaking measurement
% 
% Args:
%     As (double):                          split parameter
%     V_d_A (function):                     function of potential depending on A
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
%     N_cases (int):                        number of cases observed
%     u0_list (array):                      N_x x N_cases initial states for all cases
%     dt (double):                          time discretization
%     indx_offset (double):                 number of spacial discretization steps the ground state is offset
%
% Returns:
%     sim_freq (array):                     1  x N_cases frequency of interwell distance
%     sim_amp (array):                      1  x N_cases amplitude of interwell distance
%     sim_freq_x0 (array):                  1  x N_cases frequency of center of mass
%     sim_freq_x0_both (array):             2  x N_cases freuqencies of two highest amplitudes of center of mass sorted
%     sort_indx (int):                      index for sorting both frequencies
%     last_out_struct (struct):             contains some data of last case
%     f_fft_x0 (array):                     frequency vector
%     Y_fft_x0 (array):                     absolute fourier values 
%     sim_freq_x0_both_maxk (double):       2  x N_cases freuqencies of two highest amplitudes of center of mass unsorted
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('indx_offset','var')
    indx_offset = 2;
end

T_before_sloshing   = 0.1;
T_sloshing          = 5;

sim_freq    = zeros(1,N_cases);
sim_amp     = sim_freq;
sim_freq_x0 = sim_freq;
sim_freq_x0_both = zeros(2,N_cases);
sim_freq_x0_both_maxk = zeros(2,N_cases);


t_sloshing                                  = T_before_sloshing+dt:dt:(T_before_sloshing+T_sloshing);
t_ramp                                      = 0:dt:(T_before_sloshing);
N_t_ramp                                    = numel(t_ramp);
t                                           = [t_ramp,t_sloshing];
N_t_plot                                    = numel(t);

last_out_struct.t                           = t;



for iter_cases = 1:N_cases
    u0_shifted                                  = zeros(size(x));
    u0_shifted(indx_offset:end)                 = u0_list(1:end-indx_offset+1,iter_cases);
    A_t                                         = ones(1,N_t_plot)*As(1,iter_cases);
    psi_iter                                    = calc_prop_GPE(N_t_plot,x,u0_shifted,dt,V_d_A,A_t,par);

    [relativ_difference_iter,~,x0_iter]         = calculate_difference_from_fit(x,psi_iter);
    relative_diff_sloshing                      = relativ_difference_iter(N_t_ramp+1:end);
    x0_iter_sloshing                            = x0_iter(N_t_ramp+1:end);

    subsample_factor                            = 50;
    fft_window_factor                           = 5e3;
    relative_diff_sloshing_subsampled           = relative_diff_sloshing(1:subsample_factor:end)-mean(relative_diff_sloshing);
    N_t_sloshing_subsampled                     = numel(relative_diff_sloshing_subsampled);
    N_t_sloshing_subsampled_scaled              = 2^nextpow2(fft_window_factor*N_t_sloshing_subsampled);
    Y                                           = fft(relative_diff_sloshing_subsampled,N_t_sloshing_subsampled_scaled)/N_t_sloshing_subsampled;
    f                                           = 1/dt/subsample_factor/N_t_sloshing_subsampled_scaled*(0:N_t_sloshing_subsampled_scaled-1);
    [amp,max_indx_f]                            = max(2*abs(Y(2:round(end/2))));
    sim_freq(iter_cases)                        = f(max_indx_f+1);
    sim_amp(iter_cases)                         = amp;

    x0_sloshing_subsampled                      = x0_iter_sloshing(1:subsample_factor:end)-mean(x0_iter_sloshing(1:subsample_factor:end));
    Y                                           = fft(x0_sloshing_subsampled,N_t_sloshing_subsampled_scaled)/N_t_sloshing_subsampled;
    f                                           = 1/dt/subsample_factor/N_t_sloshing_subsampled_scaled*(0:N_t_sloshing_subsampled_scaled-1);
    % [~,max_indx_f]                              = max(2*abs(Y(2:round(end/2))));

    Y_half                                      = 2*abs(Y(2:round(end/2)));
    indices_local_max                           = islocalmax(Y_half,2);
    Y_half_diff_reduced                         = Y_half;
    Y_half_diff_reduced(~indices_local_max)     = 0;
    
    [amp_maxk,maxk_indx_f]                      = maxk(Y_half_diff_reduced,2);
    if (amp_maxk(1)-amp_maxk(2))/amp_maxk(1) < 0.4
        [max_indx_f,indx_is]                        = min(maxk_indx_f);
        [maxk_indx_f_sorted,sort_indx]              = sort(maxk_indx_f);
    else
        max_indx_f                                  = maxk_indx_f(1);
        [maxk_indx_f_sorted,sort_indx]              = sort(maxk_indx_f);
        indx_is                                     = 1;
    end
    
    sim_freq_x0(iter_cases)                     = f(max_indx_f+1);
    sim_freq_x0_both(:,iter_cases)              = f(maxk_indx_f_sorted+1);
    sim_freq_x0_both_maxk(:,iter_cases)         = f(maxk_indx_f+1);
    last_out_struct.psi_iter                    = psi_iter;
    last_out_struct.relativ_difference_iter     = relativ_difference_iter;
    last_out_struct.x0_iter                     = x0_iter;
    

    f_fft_x0                                    = f(1:end/2);
    Y_temp                                      = fft(hamming(N_t_sloshing_subsampled)'.*x0_sloshing_subsampled,N_t_sloshing_subsampled_scaled)/N_t_sloshing_subsampled;
    Y_fft_x0                                    = Y_temp(1:end/2);

end

end