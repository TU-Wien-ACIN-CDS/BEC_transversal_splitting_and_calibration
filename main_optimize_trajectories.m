addpath("functions\")
addpath("Parameters\")
addpath("Storage\Sensitivities\")

run Parameters_Transversal_Splitting

%% select calibrated model parameters

calibrated_parameter_file       = 'calibrated_parameters_March_2025.mat'; % change this to newest calibration
path_to_calibrated_experiments  = fullfile(pwd,"Storage/Calibrated_parameters",calibrated_parameter_file);

load(path_to_calibrated_experiments,'calib_struct') % loads parameters
fprintf('loaded %s\n',calibrated_parameter_file)

%% select scenario
switch_scenario = 1;    % 1: linear ramps for defined times and A_end, 
                        % 2: use times and A_end from a previous scenario, 
                        % 3: use previous scenarios ramps as initial ramps, 
                        % 4: use simulated starting point instead of ground state, A_0 values are automatically set by the selected file!!

switch_reuse    = 0;    % 0: use defined initial ramps based on switch_scenario, 
                        % 1: rescale solution of previous iteration as initial ramp for iteration

% Trehshold for energy over groundstate in normalization:
val_opt_thresh          = 1e-3; % threshold value of final energy for end of optimization, standard is 1e-3;

% IF switch_scenario == 1 or 4:
A_end_list      = [0.5];
A_0_list        = 0.28*ones(size(A_end_list));
T_list          = [1]; % in ms
N_ramps         = numel(A_end_list);

% IF switch_scenario == 2 or 3
previous_scenario_file = 'Optimized_ramps_25_11_2025_09_35_12.mat'; % file created by this code of ramps
if switch_scenario == 2 || switch_scenario == 3
    load(fullfile("Storage\Optimized_Ramps\",previous_scenario_file),"optimized_trajectories_struct")
    N_ramps = numel(optimized_trajectories_struct.ramps_cell);
    previous_struct = optimized_trajectories_struct;
else
    previous_struct = 'none';
end

% if switch_scenario == 4
ramp_to_initial_state_file = 'state_of_ramp_Nr2_after_T_hold_1.000ms_of_Optimized_ramps_with_lin_25_11_2025_09_35_12.mat'; % file create by the main_simulate_ramps code
if switch_scenario == 4
    load(fullfile(pwd,"Storage\Simulated_Initial_States\",ramp_to_initial_state_file),"initial_state_struct")
end




%% run trajectory optimization

run optimize_trajectories_internal


%% save results
optimized_trajectories_struct_without_lin.calib_struct_used                 = calib_struct;    % used calibration struct
optimized_trajectories_struct_without_lin.switch_scenario                   = switch_scenario; % selection of scenario
optimized_trajectories_struct_without_lin.switch_reuse                      = switch_reuse;    % selection for reuse variant
optimized_trajectories_struct_without_lin.previous_scenario_if_applicable   = previous_struct; % previous scenario if applicable
optimized_trajectories_struct_without_lin.ramps_cell                        = {};
for indx_ramp = 1:N_ramps
    ramp_struct.is_ramp = 1;
    if size(A_t_is_cases{indx_ramp},1)>size(A_t_is_cases{indx_ramp},2)
        ramp_struct.A_t     = A_t_is_cases{indx_ramp}';
    else
        ramp_struct.A_t     = A_t_is_cases{indx_ramp};
    end
    ramp_struct.t       = t_case{indx_ramp};
    ramp_struct.A_0     = ramp_struct.A_t(1);
    ramp_struct.A_end   = ramp_struct.A_t(end);
    ramp_struct.T_ramp  = ramp_struct.t(end);
    ramp_struct.A_initi = A_t_init_used{indx_ramp};
    optimized_trajectories_struct_without_lin.ramps_cell{indx_ramp} = ramp_struct;
end
optimized_trajectories_struct       = optimized_trajectories_struct_without_lin; % ramps without linear ramps inbetween

timestampStr                        = char(datetime('now', 'Format', 'dd_MM_yyyy_HH_mm_ss'));
uniqueFileName_optimized_ramps      = ['Optimized_ramps_', timestampStr, '.mat'];
uniqueFilePath_optimized_ramps      = fullfile(pwd,"Storage/Optimized_Ramps/", uniqueFileName_optimized_ramps);
save(uniqueFilePath_optimized_ramps,'optimized_trajectories_struct')
fprintf("ramps without linear saved to: %s\n",fullfile("Storage/Optimized_Ramps/", uniqueFileName_optimized_ramps))


optimized_trajectories_struct_with_lin.calib_struct_used                = calib_struct;    % used calibration struct
optimized_trajectories_struct_with_lin.switch_scenario                  = switch_scenario; % selection of scenario
optimized_trajectories_struct_with_lin.switch_reuse                     = switch_reuse;    % selection for reuse variant
optimized_trajectories_struct_with_lin.previous_scenario_if_applicable  = previous_struct; % previous scenario if applicable
optimized_trajectories_struct_with_lin.ramps_cell                       = {};
for indx_ramp = 1:N_ramps
    ramp_struct     = optimized_trajectories_struct_without_lin.ramps_cell{indx_ramp};
    ramp_struct.A_t = linspace(ramp_struct.A_0,ramp_struct.A_end,numel(ramp_struct.A_t));
    optimized_trajectories_struct_with_lin.ramps_cell{2*(indx_ramp-1)+1} = optimized_trajectories_struct_without_lin.ramps_cell{indx_ramp};
    optimized_trajectories_struct_with_lin.ramps_cell{2*(indx_ramp-1)+2} = ramp_struct;
end
optimized_trajectories_struct       = optimized_trajectories_struct_with_lin; % ramps with linear ramps inbetween

uniqueFileName_optimized_ramps      = ['Optimized_ramps_with_lin_', timestampStr, '.mat'];
uniqueFilePath_optimized_ramps      = fullfile(pwd,"Storage/Optimized_Ramps/", uniqueFileName_optimized_ramps);
save(uniqueFilePath_optimized_ramps,'optimized_trajectories_struct')
fprintf("ramps with linear saved to: %s\n",fullfile("Storage/Optimized_Ramps/", uniqueFileName_optimized_ramps))

% if starting points are from previous simulation
if switch_scenario == 4
    optimized_trajectories_struct_sim_initial.calib_struct_used                 = calib_struct;    % used calibration struct
    optimized_trajectories_struct_sim_initial.switch_scenario                   = switch_scenario; % selection of scenario
    optimized_trajectories_struct_sim_initial.switch_reuse                      = switch_reuse;    % selection for reuse variant
    optimized_trajectories_struct_sim_initial.previous_scenario_if_applicable   = previous_struct; % previous scenario if applicable
    optimized_trajectories_struct_sim_initial.ramps_cell                        = {};
    optimized_trajectories_struct_sim_initial.initial_state_struct              = initial_state_struct; % struct containing previous ramp and initial state
    for indx_ramp = 1:N_ramps
        ramp_struct.is_ramp = 1;
        if size(A_t_is_cases{indx_ramp},1)>size(A_t_is_cases{indx_ramp},2)
            ramp_struct.A_t     = A_t_is_cases{indx_ramp}';
        else
            ramp_struct.A_t     = A_t_is_cases{indx_ramp};
        end
        ramp_struct.t       = t_case{indx_ramp};
        ramp_struct.A_0     = ramp_struct.A_t(1);
        ramp_struct.A_end   = ramp_struct.A_t(end);
        ramp_struct.T_ramp  = ramp_struct.t(end);
        ramp_struct.A_initi = A_t_init_used{indx_ramp};
        optimized_trajectories_struct_sim_initial.ramps_cell{indx_ramp} = ramp_struct;
    end
    optimized_trajectories_struct       = optimized_trajectories_struct_sim_initial;
    uniqueFileName_optimized_ramps      = ['Optimized_ramps_from_initial_sim_', timestampStr, '.mat'];
    uniqueFilePath_optimized_ramps      = fullfile(pwd,"Storage/Optimized_Ramps/", uniqueFileName_optimized_ramps);
    save(uniqueFilePath_optimized_ramps,'optimized_trajectories_struct')
    fprintf("ramps with initial data saved to: %s\n",fullfile("Storage/Optimized_Ramps/", uniqueFileName_optimized_ramps))

end
