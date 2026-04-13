function [J]    = cost_ramps_and_offset_only_f(params,V_fun,x,par_orig,dt,exp_ramps,exp_offset,u0_ramp_init,u0_offset_init,f_ramps_exp,f_offset_exp,par_norm,offset_shift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates cost from ramp selection and measurements
% 
% Args:
%     params (double):                      split parameter
%     V_fun (function):                     potential from paramters, A and x
%     x (array):                            1 x N_x spatial values
%     par_orig (struct):                    model parameters before adaptation from BFGS
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
%     dt (double):                          time discretization
%     exp_ramps (struct):                   ramps from calibration_experiment_struct
%     exp_offset (struct):                  offset measurements from calibration_experiment_struct
%     u0_ramp_init (array):                 N_x x N_ramps ground states from initial parameters as initial guesses for ground states
%     u0_offset_init (array):               N_x x N_offset ground states from initial parameters as initial guesses for ground states
%     f_ramps_exp (array):                  1/2 x N_ramps experimental measurements of ramps used for calibration
%     f_offset_exp (array):                 1/2 x N_offset experimental measurements of offset measurements used for calibration
%     par_norm (array):                     1 x N_parameters normalization values for parameters of potential and g1D/g_perp
%     offset_shift (int array):             1 x N_offset number of spacial discretization steps the ground state is offset
%
% Returns:
%     J (double):                           cost function value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_sim       = @(A)V_fun(params,A,x);

par     = par_orig;
par.g   = params(end)*par_norm(end);


u0p28  = calculate_groundstate_imag_time_GPE_ss(V_sim(0.28),par,u0_ramp_init,x);

N_offset    = size(exp_offset,2);
N_ramps     = size(exp_ramps,2);

f_ramps_sim    = zeros(size(f_ramps_exp));

f_offset_init  = zeros(size(f_offset_exp));

u0_offset   = zeros(numel(x),N_offset);
for indx_offset     = 1:N_offset
    u0_offset(:,indx_offset)    = calculate_groundstate_imag_time_GPE_ss(V_sim(exp_offset{indx_offset}.A_offset),par,u0_offset_init(:,indx_offset),x);
end

for indx_ramp = 1:N_ramps
    sloshy                      = freq_and_amp_from_ramp(exp_ramps{indx_ramp}.A_t,dt,V_sim,x,par,u0p28,3);
    f_ramps_sim(1,indx_ramp)    = sloshy(1);

    if f_ramps_exp(1,indx_ramp) ~= f_ramps_exp(2,indx_ramp)
        f_ramps_sim(:,indx_ramp)    = sort(sloshy([1,11]));
    else
        f_ramps_sim(:,indx_ramp)    = sloshy(1);
    end
end

for indx_offset = 1:N_offset
    [~,~,x0_freq,x0_freq_both,sort_indx]    = calculate_shaking_frequencies(exp_offset{indx_offset}.A_offset,V_sim,x,par,1,u0_offset(:,indx_offset),dt,offset_shift(indx_offset));
    if size(f_offset_init,1) == 1
        f_offset_init(:,indx_offset)        = x0_freq;
    else
        if f_offset_exp(1,indx_offset) == f_offset_exp(2,indx_offset)
            f_offset_init(:,indx_offset)    = x0_freq_both(sort_indx == 1);
        else
            f_offset_init(:,indx_offset)    = x0_freq_both;
        end
    end
end


J_ramps     = sum((f_ramps_exp - f_ramps_sim).^2);
J_offset    = sum((f_offset_init - f_offset_exp).^2);

J   = J_ramps + J_offset;

% fprintf("%.9e\n",params')
fprintf("%.6e + %.6e = %.6e\n",J_ramps,J_offset,J)

end