addpath("functions\")
addpath("Parameters\")
addpath("Storage\Sensitivities\")

run Parameters_Transversal_Splitting


%% estimated model parameters
k1_est                  = -1495;             % slope of curvature before As
k2_est                  = 2210;              % slope of curvature after As
c_est                   = 2.115;                % factor before sqrt function for minimum position c*sqrt(A-As)
As_est                  = 0.4;                  % parameter at which the double well starts forming
g1D_est                 = g1D_first_principle;  % nonlinear interaction constant for the transversal direction


%% IMPORTANT: selection if GA or heuristic method is used: 
% GA needs to precalculate the sensitivities, can take a long time (will reuse data if possible) but optimizes the information content of the experiment selection, 
% heuristic is not optimal but expected to be relatively good, based on patterns in results of GA, only depends on value of As_est and range of possible experiments
use_heuristic_method    = 1;    % 1: heuristic method is used, 0: GA is used,


%% parameters for the calibration experiment selection

% range of possible experiments
A_ramp_min      = 0.50;         % minimal value at which double well is visible
A_ramp_max      = 0.65;         % maximal value where good measurements are possible
A_offset_min    = 0.28;         % minimal value where offset measurements are used (generally value used for initial value of splitting)
T_ramp_min      = 1;            % minimal time for linear ramps (must produce measuremable fringe for A_ramp_max)
T_ramp_max      = 3;            % maximal time for linear ramps
A_ramps_start   = A_offset_min; % starting point of ramps


% hyperparameter for the genetic algorithm
CULLING_Variant         = 1;    % 1: cull the 2/3 of population with highest cost, 2: cull random 2/3 of population, 3: use linear propability distribution for culling
mutation_var            = 1;    % 1: mutation of all members, 2: mutation of only new members


%%

% run the desired parameters, should not be tampered with!
run experiment_selection_internal


%% Save results uniquely identifiable
% get unique timestamp
timestampStr                = char(datetime('now', 'Format', 'dd_MM_yyyy_HH_mm_ss'));
uniqueFileName_ramps_offset = ['ramps_and_offset_for_Experiment_', timestampStr, '.mat'];
uniqueFileName_selection    = ['selection_', timestampStr, '.mat'];

% save results
uniqueFilePath_ramps_offset = fullfile(pwd,"Storage/Calibration_Szenarios/", uniqueFileName_ramps_offset);
save(uniqueFilePath_ramps_offset,'calibration_experiment_struct')
uniqueFilePath_selection    = fullfile(pwd,"Storage/Calibration_Szenarios/", uniqueFileName_selection);
save(uniqueFilePath_selection,'selection_struct')

% display paths to results
disp('Selection and ramps/offset values are saved as')
disp(uniqueFilePath_ramps_offset)
disp(uniqueFilePath_selection)