addpath("functions\")
addpath("Parameters\")
addpath("Storage\Sensitivities\")
% set(0,'DefaultFigureWindowStyle','docked')

run Parameters_Transversal_Splitting

%% select calibrated model parameters

calibrated_parameter_file       = 'calibrated_parameters_March_2025.mat'; % change this to newest calibration
path_to_calibrated_experiments  = fullfile(pwd,"Storage/Calibrated_parameters",calibrated_parameter_file);

load(path_to_calibrated_experiments,'calib_struct') % loads parameters
fprintf('loaded %s\n',calibrated_parameter_file)

%% select scenario
optimized_trajectories_file     = 'Optimized_ramps_with_lin_18_03_2026_10_04_43.mat'; % change this to optimized ramps of scenario
path_to_optimized_trajectories  = fullfile(pwd,'Storage/Optimized_Ramps',optimized_trajectories_file);

load(path_to_optimized_trajectories,'optimized_trajectories_struct')

%% select time for waveform simulation of certain ramp for use as initial state of further optimization
desired_hold_Time       = 1; % in ms
ramp_selected           = 2; % if ramp_selected > N_ramp -> ramp selected will be set to N_ramp

%% run scenarios to get results

run simulate_trajectory_internal


%% save results

results_file_name               = strcat('results_of_',optimized_trajectories_file);
path_to_results_file            = fullfile(pwd,"Storage/Simulated_Data",results_file_name);
save(path_to_results_file,"GPE_results_struct","-v7.3")
fprintf('results saved to %s\n',results_file_name)

%% save result of ramp at desired point
initial_state_struct_file_name          = strcat(sprintf('state_of_ramp_Nr%i_after_T_hold_%.3fms_of_',ramp_selected,desired_hold_Time),optimized_trajectories_file);
path_to_initial_state_struct_file       = fullfile(pwd,"Storage/Simulated_Initial_States",initial_state_struct_file_name);
save(path_to_initial_state_struct_file,"initial_state_struct","-v7.3")
fprintf('initial_state_struct saved to %s\n',initial_state_struct_file_name)

%% plot results
plot_experiment_data_results(GPE_results_struct)