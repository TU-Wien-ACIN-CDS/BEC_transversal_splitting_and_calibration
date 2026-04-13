% add path to functions and parameters
addpath("functions\")
addpath("Parameters\")

% loading parameters
run Parameters_Transversal_Splitting

%% Loading structure from experiment selection, as well as measured frequencies

% struct with ramps and offset values
calibration_experiment_struct_file_name     = 'ramps_and_offset_for_Experiment_25_11_2025_09_26_37.mat';
path_to_calibration_experiment_struct       = fullfile(pwd,'Storage','Calibration_Szenarios',calibration_experiment_struct_file_name);
load(path_to_calibration_experiment_struct,'calibration_experiment_struct')

% measurement values, make sure the sequence of measurements is the same as in the struct with ramps and offset values above
f_values_ramps_exp  = [;];  % add in values of experiment measurements, can be two rows, 
                            % if column entries are different, two frequencies are use for THIS ramp (not all), put higher frequency in second row!
                            % [1,1.5,2;1,3.5,2] in kHz
f_values_offset_exp = [;];  % add in values of experiment measurements, can be two rows, 
                            % if column entries are different, two frequencies are use for THIS value of A (not all), put higher frequency in second row!



%% Parameters for calibration
use_estimates_from_selection = 1;   % 1: initial estimates for model parameters are taken from values used for selection of experiments, 0: original fitted parameters are used. 
                                    % select 0 if "calibration_experiment_struct" does not come from the "main_experiment_selection" code!



run calibrate_model_internal

%%

% creating unique timestamp
timestampStr                        = char(datetime('now', 'Format', 'dd_MM_yyyy_HH_mm_ss'));
uniqueFileName_calib_struct         = ['calibrated_parameters_', timestampStr, '.mat'];

% storing calibrated parameters
calib_struct.k1                     = k1_opt;  % slope of curvature before As
calib_struct.k2                     = k2_opt;  % slope of curvature after As
calib_struct.c                      = c_opt;   % factor before sqrt function for minimum position x_m = c*sqrt(A-As)
calib_struct.As                     = As_opt;  % parameter at which the double well starts forming
calib_struct.g1D                    = g1D_opt; % nonlinear interaction constant for the transversal direction
calib_struct.coeffs_a6              = coeffs_a6_opt; % coefficients of a6(A)*x^6 polynomial
calib_struct.coeffs_a4              = coeffs_a4_opt; % coefficients of a4(A)*x^4 polynomial
calib_struct.coeffs_a2_pre_As       = coeffs_a2_pre_As_opt;  % coefficients of a2(A)*x^2 polynomial before As
calib_struct.coeffs_a2_post_As      = coeffs_a2_post_As_opt; % coefficients of a2(A)*x^2 polynomial after As
calib_struct.coeffs_a0              = coeffs_a0_opt; % coeffiecients of constant part of polynomial
calib_struct.f_values_ramps         = f_values_ramps_exp;  % used frequency measurements of ramps in calibration
calib_struct.f_values_offset        = f_values_offset_exp; % used frequency measurements of offset in calibration
calib_struct.timestamp              = timestampStr; % timestamp used

uniqueFilePath_calib_struct         = fullfile(pwd,"Storage/Calibrated_parameters/", uniqueFileName_calib_struct);
save(uniqueFilePath_calib_struct,'calib_struct')
fprintf('calibrated_parameters where saved under: %s\n',uniqueFilePath_calib_struct)