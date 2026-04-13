function [flag_recompute,location,paramFilePath] = check_if_sensitivities_already_exist(k1_est,k2_est,c_est,As_est,g1D_est,A_ramp_min,A_ramp_max,T_ramp_min,T_ramp_max,A_offset_min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if sensitivity needs to be recalculated (partially written by Gemini (LLM))
% 
% Args:
%     k1_est (double):                      estimated slope of curvature before split parameter
%     k2_est (double):                      estimated slope of curvature after split parameter
%     c_est (double):                       estimated factor before sqrt function for minimum position x_m = c*sqrt(A-As)
%     As_est (double):                      estimated split parameter
%     g1D_est (double):                     estimated g_perp
%     A_ramp_min (double):                  minimum ramp value allowed
%     A_ramp_max (double):                  maximum ramp value allowed
%     T_ramp_min (double):                  minimum ramp time allowed
%     T_ramp_max (double):                  maximum ramp time allowed
%     A_offset_min (double):                minimum value for offset measurements
%
% Returns:
%     flag_recompute (int):                 if 1 sensitivities need to be recomputed
%     location (string):                    path to sensitivities
%     paramFilePath (string):               path to paramter file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


parameter_change_threshold = 0.01;
path_to_storage = "Storage/Sensitivities";

parameterFileName = 'used_parameters.mat';

searchPattern = 'V\d+';

currentDirContents = dir(path_to_storage);

dirNames = {currentDirContents.name};
isDir = [currentDirContents.isdir];

versionDirs = {};
versionNumbers = [];

for i = 1:length(dirNames)
    if isDir(i) && ~isempty(regexp(dirNames{i}, ['^' searchPattern '$'], 'once'))
        versionDirs{end+1} = dirNames{i};
        % Extract the numeric part (e.g., '1' from 'V1')
        versionNumbers(end+1) = str2double(regexprep(dirNames{i}, 'V', ''));
    end
end

% --- 2. Iterate through existing folders and check parameters ---
% Sort folders numerically to ensure V1 is checked before V10
[~, sortIdx] = sort(versionNumbers);

for i = 1:length(sortIdx)
    currentFolderName = versionDirs{sortIdx(i)};
    currentFolderPath = fullfile(pwd,path_to_storage, currentFolderName);
    paramFilePath = fullfile(pwd,path_to_storage,currentFolderName, parameterFileName);
    
    fprintf('Checking folder: %s...\n', currentFolderName);

    % Check if the parameter file exists
    if exist(paramFilePath, 'file') == 2
        try
            % Load the parameter file content
            load(paramFilePath,'data');
            
            % Check if the comparison field exists and is numeric

            check_k1    = ((data.k1-k1_est)/data.k1)<parameter_change_threshold;
            check_k2    = ((data.k2-k2_est)/data.k2)<parameter_change_threshold;
            check_c     = ((data.c-c_est)/data.c)<parameter_change_threshold;
            check_As    = ((data.As-As_est)/data.As)<parameter_change_threshold;
            check_g1D   = ((data.g1D-g1D_est)/data.g1D)<parameter_change_threshold;
            % Perform the comparison using tolerance (crucial for floats)
            if check_k1 && check_k2 && check_c && check_As && check_g1D && (data.A_offset_min == A_offset_min) && (data.A_ramp_max == A_ramp_max) && (data.A_ramp_min == A_ramp_min) && (data.T_ramp_max == T_ramp_max) && (data.T_ramp_min == T_ramp_min)
                % MATCH FOUND!
                if ~data.is_done
                    disp('Data at')
                    disp(currentFolderPath)
                    disp('is not fully calculated!')
                    break
                end
                disp('>> MATCH FOUND!');
                location        = currentFolderPath;
                flag_recompute  = 0;
                return; % Exit the function and return the matching path
            end
        catch ME
            fprintf('   Error loading or accessing file %s: %s\n', parameterFileName, ME.message);
        end
    else
        fprintf('   File %s not found in this folder.\n', parameterFileName);
    end
end   

% --- 3. No match found: Create a new version folder ---
disp('--- No matching version found. Creating a new folder ---');

if isempty(versionNumbers)
    nextVersion = 1; % Start at V1 if no V folders exist
else
    nextVersion = max(versionNumbers) + 1; % Increment max version number
end

newFolderName = ['V', num2str(nextVersion)];
newFolderPath = fullfile(pwd,path_to_storage, newFolderName);

% Create the directory
[success, msg, msgID] = mkdir(newFolderPath);

data.k1             = k1_est;
data.k2             = k2_est;
data.c              = c_est;
data.As             = As_est;
data.g1D            = g1D_est;
data.A_ramp_min     = A_ramp_min;
data.A_ramp_max     = A_ramp_max;
data.T_ramp_min     = T_ramp_min;
data.T_ramp_max     = T_ramp_max;
data.A_offset_min   = A_offset_min;
data.is_done        = false;

paramFilePath   = fullfile(newFolderPath,parameterFileName);
save(paramFilePath,'data')


if success
    fprintf('Successfully created new folder: %s\n', newFolderName);
    folderPath = newFolderPath;
else
    warning('Failed to create new folder %s: %s (%s)', newFolderName, msg, msgID);
    folderPath = ''; % Return empty on failure
end

flag_recompute  = 1;
location        = folderPath;
end