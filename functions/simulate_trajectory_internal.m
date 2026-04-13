%% This is an internal function, be very carefull when changing anything!!!!

%% simulation parameters
T_sloshing          = 3; % length of sloshing phase in ms
k_treshold          = 25; % range of k values considered in momentum calculated

if T_sloshing < desired_hold_Time
    T_sloshing = desired_hold_Time;
end

%% Defining potential
parameters_calibrated   = [calib_struct.k1,calib_struct.k2,calib_struct.c,calib_struct.As]; % calibrated parameters from calibration file
par.g                   = calib_struct.g1D; % calibrated g1D
g1D                     = par.g;
m                       = par.m;

x_bound             = 2.5;  % max range of potential
a6_coeffs           = 0;    % coefficients of a6(A)*x6 part
a4_coeffs           = @(x_is)x_is(2)/8/x_is(3)^2; % coefficients of a4(A)*x4 part
a2_coeffs_pre_As    = @(x_is)[x_is(1)/2, -x_is(1)/2*x_is(4)];  % coefficients of a2(A)*x2 part before As
a2_coeffs_post_As   = @(x_is)[- x_is(2)/4, x_is(2)/4*x_is(4)]; % coefficients of a2(A)*x2 part after As
a2_fun              = @(x_is,A) polyval(heaviside(x_is(4)-A)*a2_coeffs_pre_As(x_is) + heaviside(A-x_is(4))*a2_coeffs_post_As(x_is),A);  % combined function of a2(A)*x^2
V_fun               = @(x_is,A,x) polyval([a4_coeffs(x_is),0,a2_fun(x_is,A),0,0],x).*heaviside(x+x_bound).*heaviside(x_bound-x)+1e5.*heaviside(-x-x_bound)+1e5.*heaviside(-x_bound+x); % combined V as function x and parameters
V_A                 = @(A) V_fun(parameters_calibrated,A,x); % Potential as function of A

ham                 = @(A)-lap4/(2*m)+spdiag(V_A(A)); % hamiltonian component
energy_fun          = @(u,A) real(trapz(x,conj(u).*(ham(A)+par.g/2*(spdiag(abs(u).^2)))*u)); % function to calculate energy



%% Calculate initial states
fprintf('Calculating initial and final ground states...\n')
N_ramps     = numel(optimized_trajectories_struct.ramps_cell);

if ramp_selected > N_ramps
    ramp_selected = N_ramps;
end

u0_case     = zeros(N_x,N_ramps);
u1_case     = zeros(N_x,N_ramps);
for indx_ramp = 1:N_ramps
    if indx_ramp == 1 % first ground state
        u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(optimized_trajectories_struct.ramps_cell{indx_ramp}.A_0),par,ones(size(x)),x);
        u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(optimized_trajectories_struct.ramps_cell{indx_ramp}.A_end),par,ones(size(x)),x);
    else % after first ground state
        if optimized_trajectories_struct.ramps_cell{indx_ramp}.A_0 == optimized_trajectories_struct.ramps_cell{indx_ramp-1}.A_0 % if A_0 is same as previous
            u0_case(:,indx_ramp)        = u0_case(:,indx_ramp-1);
        else
            u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(optimized_trajectories_struct.ramps_cell{indx_ramp}.A_0),par,u0_case(:,indx_ramp-1),x);
        end
        if optimized_trajectories_struct.ramps_cell{indx_ramp}.A_end == optimized_trajectories_struct.ramps_cell{indx_ramp-1}.A_end % if A_end is same as previous
            u1_case(:,indx_ramp)        = u1_case(:,indx_ramp-1);
        else
            u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(optimized_trajectories_struct.ramps_cell{indx_ramp}.A_end),par,u1_case(:,indx_ramp-1),x);
        end
    end
end
fprintf('Done\n')

%% initialize data
psi_iter                = cell(N_ramps,1);
t_case                  = cell(N_ramps,1);
A_t_case                = cell(N_ramps,1);
trap_bottom_case        = cell(N_ramps,1);
energy_t                = cell(N_ramps,1);
rel_diff_sim            = cell(N_ramps,1);
k_sim                   = cell(N_ramps,1);
m_sim                   = cell(N_ramps,1);
t_sim_momentum          = cell(N_ramps,1);
sigma_sim               = cell(N_ramps,1);
A_t_moment              = cell(N_ramps,1);

f_ramps_sim             = zeros(N_ramps,1);
amp_ramps_sim           = zeros(N_ramps,1);
d_ramps_sim             = zeros(N_ramps,1);
ph_ramps_sim            = zeros(N_ramps,1);
var_d_ramps_sim         = zeros(N_ramps,1);
mean_sigma_sim          = zeros(N_ramps,1);
var_sigma_ramps_sim     = zeros(N_ramps,1);
amp_sigma_sim           = zeros(N_ramps,1);
f_sigma_ramps_sim       = zeros(N_ramps,1);
var_sigma_M_ramps_sim   = zeros(N_ramps,1);
T_ramp                  = zeros(N_ramps,1);
A_end_ramp              = zeros(N_ramps,1);
energy_sim              = zeros(N_ramps,1);
energy_GS_end           = zeros(N_ramps,1);

%% running simulation
fprintf('Starting simulation of ramps...\n')

for indx_ramp = 1:N_ramps
    % get time discretization
    dt                                  = optimized_trajectories_struct.ramps_cell{indx_ramp}.t(2) - optimized_trajectories_struct.ramps_cell{indx_ramp}.t(1);
    % calculate simulation result of trajectory
    [sloshy,psi_iter{indx_ramp},A_t_case{indx_ramp},t_case{indx_ramp}]    = freq_and_amp_from_ramp(optimized_trajectories_struct.ramps_cell{indx_ramp}.A_t,dt,V_A,x,par,u0_case(:,indx_ramp),T_sloshing);

    % extract data from sloshy
    f_ramps_sim(indx_ramp)              = sloshy(1); % frequency of interwell distance
    amp_ramps_sim(indx_ramp)            = sloshy(2); % amplitude of interwell distance
    d_ramps_sim(indx_ramp)              = sloshy(3); % mean of interwell distance
    ph_ramps_sim(indx_ramp)             = sloshy(4); % phase of interwell distance
    var_d_ramps_sim(indx_ramp)          = sloshy(5); % variane of interwell distance
    mean_sigma_sim(indx_ramp)           = sloshy(6); % mean of sigma in situ
    var_sigma_ramps_sim(indx_ramp)      = sloshy(7); % variance of sigma in situ
    amp_sigma_sim(indx_ramp)            = sloshy(8); % amplitude of sigma in situ
    f_sigma_ramps_sim(indx_ramp)        = sloshy(9); % frequency of sigma in situ
    var_sigma_M_ramps_sim(indx_ramp)    = sloshy(10); % variance of approximated sigma in momentum space (using sigma_M = 1/(sqrt(2)*sigam_Psi) )
    T_ramp(indx_ramp)                   = optimized_trajectories_struct.ramps_cell{indx_ramp}.t(end); % ramp length
    A_end_ramp(indx_ramp)               = optimized_trajectories_struct.ramps_cell{indx_ramp}.A_end;  % end value of trajectory
    energy_sim(indx_ramp)               = energy_fun(psi_iter{indx_ramp}(:,end),A_end_ramp(indx_ramp)); % energy at end of simulation
    energy_GS_end(indx_ramp)            = energy_fun(u1_case(:,indx_ramp),A_end_ramp(indx_ramp));       % energy of ground state of A_end

    [rel_diff_sim{indx_ramp},sigma_sim{indx_ramp}]         = calculate_difference_from_fit(x,psi_iter{indx_ramp}); % relative difference for current trajectory

    [k_sim{indx_ramp},m_sim{indx_ramp}] = calc_moment(x,psi_iter{indx_ramp},16*N_x,1,k_treshold); % momentum representation for current trajectory

    fprintf('(%i/%i) done\n',indx_ramp,N_ramps)
end


%% collecting result struct

GPE_results_struct.x                    = x;
GPE_results_struct.calib_struct_used    = calib_struct;
GPE_results_struct.N_ramps              = N_ramps;
GPE_results_struct.Szenario_file_name   = optimized_trajectories_file;
Szenarios = cell(N_ramps,1);
for indx_ramp = 1:N_ramps
    Szenario.t              = t_case{indx_ramp};
    Szenario.A_t            = A_t_case{indx_ramp};
    Szenario.T_ramp         = T_ramp(indx_ramp);
    Szenario.Psi_x_t        = psi_iter{indx_ramp};
    Szenario.d_t            = rel_diff_sim{indx_ramp};
    Szenario.sigma_t        = sigma_sim{indx_ramp};
    Szenario.d_amp          = amp_ramps_sim(indx_ramp);
    Szenario.sigma_amp      = amp_sigma_sim(indx_ramp);
    Szenario.d_freq         = f_ramps_sim(indx_ramp);
    Szenario.sigma_freq     = f_sigma_ramps_sim(indx_ramp);
    Szenario.d_mean         = d_ramps_sim(indx_ramp);
    Szenario.sigma_mean     = mean_sigma_sim(indx_ramp);
    Szenario.d_var          = var_d_ramps_sim(indx_ramp);
    Szenario.sigma_var      = var_sigma_ramps_sim(indx_ramp);
    Szenario.delta_E        = energy_sim(indx_ramp) - energy_GS_end(indx_ramp);
    Szenario.t_moment       = t_sim_momentum{indx_ramp};
    Szenario.k              = k_sim{indx_ramp};
    Szenario.Psi_hat_k_t    = m_sim{indx_ramp};
    Szenarios{indx_ramp}    = Szenario;
end
GPE_results_struct.Szenarios = Szenarios;


%% creating initial state struct for further optimization

indx_time_selected  = find(Szenarios{ramp_selected}.t > optimized_trajectories_struct.ramps_cell{ramp_selected}.T_ramp+desired_hold_Time, 1);

initial_state_struct.calib_struct_used      = calib_struct;
initial_state_struct.t                      = Szenarios{ramp_selected}.t(1:indx_time_selected);
initial_state_struct.A_t_previous           = Szenarios{ramp_selected}.A_t(1:indx_time_selected);
initial_state_struct.Psi_end                = Szenarios{ramp_selected}.Psi_x_t(:,indx_time_selected);
initial_state_struct.A_initial              = Szenarios{ramp_selected}.A_t(indx_time_selected);

