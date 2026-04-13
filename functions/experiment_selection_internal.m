%% This is an internal function, be very carefull when changing anything!!!!

A_possible_linear_ramps     = A_ramp_min:0.01:A_ramp_max; % range of A values possible for linear ramps
T_possible_linear_ramps     = T_ramp_min:0.1:T_ramp_max;  % range of T values possible for linear ramps

A_possible_CLR              = A_ramp_min+0.01:0.01:A_ramp_max; % range of A values possible for consecutive linear ramps

A_possible_Offset           = A_offset_min:0.005:A_ramp_max;   % range of A values possible for offset measurements

if use_heuristic_method
    % linear ramps: As+0.07, T = linspace(T_min,T_max,5)
    % concecutive linear ramps: As+0.07
    % Offset parameter values: [A_min,A_min+0.01,As-0.02,As-0.01,As,As+0.01,As+0.02,A_max-0.01,A_max]
    A_linear_ramps              = (As_est+0.11)*ones(1,5);
    T_linear_ramps              = linspace(T_ramp_min,T_ramp_max,5);
    
    fprintf('starting heuristic selection\n')
    A_consecutive_linear_ramps  = As_est+0.11;

    A_offset                    = [A_offset_min,A_offset_min+0.01,As_est-0.02,As_est-0.01,As_est,As_est+0.01,As_est+0.02,A_ramp_max-0.01,A_ramp_max];

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
    data.A_ramps_start  = A_ramps_start;
    data.is_done        = 'heuristic';

else
    %% Starting GA
    %% initializing all discretized spaces for GA
    fprintf('starting genetic algorithm\n')
    N_par                       = 5;                        % number of parameters to be calibrated

    A_0                         = A_ramps_start;            % initial parameter
    A_end_list                  = A_possible_linear_ramps;  % end parameters
    T_list                      = T_possible_linear_ramps;  % ramp times
    T_list                      = [T_list];                 % -||- for adding additional discrete times
    T_list                      = sort(T_list);             % sorting
    N_A                         = numel(A_end_list);        % number of end parameters
    N_T                         = numel(T_list);            % number of ramp times
    N_all                       = N_A*N_T;                  % absolute number of single linear ramps
    
    % parameters for offset sloshing measurements (O)
    A_offset_list               = A_possible_Offset;        % parameter list for offset sloshing measurements
    N_A_offset                  = numel(A_offset_list);     % number of parameters for offset sloshing measurements
    
    % parameters for multiple linear ramp experiments (ML)
    A_end_list_ML               = A_possible_CLR;               % end parameters
    A_1_list_ML                 = A_ramp_min*ones(size(A_end_list_ML)); % intermediate parameter list
    T2_list_ML                  = 5;                            % ramp time of second ramps
    T1_list_ML                  = 1*ones(size(T2_list_ML));     % ramp time of first ramps
    T_end_list_ML               = T2_list_ML + T1_list_ML;      % combined ramp times
    N_A_ML                      = numel(A_end_list_ML);         % number of end parameters
    N_T_ML                      = numel(T1_list_ML);            % number of ramp times
    N_ML                        = N_T_ML * N_A_ML;              % absolute number of multiple linear ramp experiments
    
    slosh_struct_cell           = cell(N_A,N_T);                % cells to save slosh structs
    sensitivity_f_all           = zeros(N_A,N_T,N_par);         % sensitivities of single linear ramp frequencies
    sensitivity_amp_all         = zeros(N_A,N_T,N_par);         % sensitivities of single linear ramp amplitudes
    sensitivity_d_all           = zeros(N_A,N_T,N_par);         % sensitivities of single linear ramp mean relative difference
    sensitivity_ph_all          = zeros(N_A,N_T,N_par);         % sensitivities of single linear ramp phase
    sensitivity_offset_all      = zeros(N_A_offset,N_par);      % sensitivities of offset sloshing frequencies
    sensitivity_f_ML_all        = zeros(N_A_ML,N_T_ML,N_par);   % sensitivities of multiple linear ramp frequencies
    sensitivity_amp_ML_all      = zeros(N_A_ML,N_T_ML,N_par);   % sensitivities of multiple linear ramp amplitudes
    sensitivity_ph_ML_all       = zeros(N_A_ML,N_T_ML,N_par);   % sensitivities of multiple linear ramp phase
    f_grid_all                  = zeros(N_A,N_T);               % single linear ramp frequencies
    amp_grid_all                = zeros(N_A,N_T);               % single linear ramp amplitudes
    A_end_grid_all              = zeros(N_A,N_T);               % single linear ramp end parameters
    T_grid_all                  = zeros(N_A,N_T);               % single linear ramp ramp times
    f_offset                    = zeros(N_A_offset,1);          % offset sloshing frequencies
    f_ML                        = zeros(N_A_ML,N_T_ML);         % multiple linear ramp frequencies
    amp_ML                      = zeros(N_A_ML,N_T_ML);         % multiple linear ramp amplitudes
    ph_ML                       = zeros(N_A_ML,N_T_ML);         % multiple linear ramp phases
    A_end_grid_ML               = zeros(N_A_ML,N_T_ML);         % multiple linear ramp end parameters
    T_grid_ML                   = zeros(N_A_ML,N_T_ML);         % multiple linear ramp ramp times

    %% check if sensitivities are allready calculated:
    [flag_recompute,location,path_to_parameter] = check_if_sensitivities_already_exist(k1_est,k2_est,c_est,As_est,g1D_est,A_ramp_min,A_ramp_max,T_ramp_min,T_ramp_max,A_offset_min);
    load(path_to_parameter,'data')
    % [flag_recompute,location] = check_if_sensitivities_already_exist(1,2,3,4,5,6,7,8,9,10)
    if ~flag_recompute
        fprintf('loading data\n')
        % loading single linear ramp sensitivities and values
        for indx_A = 1:N_A
            for  indx_T  = 1:N_T
                load_string                             = fullfile(location,sprintf('Sensitivity_LR_A_end_%.3f_T_%.3f.mat',A_end_list(indx_A),T_list(indx_T)));
                load(load_string,"slosh_struct")
                % disp(load_string)
                slosh_struct_cell{indx_A,indx_T}        = slosh_struct;
                S_sloshy                                = slosh_struct.S_sloshy(:,:);
                sensitivity_f_all(indx_A,indx_T,:)      = S_sloshy(1,:);
                sensitivity_amp_all(indx_A,indx_T,:)    = S_sloshy(2,:);
                sensitivity_d_all(indx_A,indx_T,:)      = S_sloshy(3,:);
                sensitivity_ph_all(indx_A,indx_T,:)     = S_sloshy(4,:);
                f_grid_all(indx_A,indx_T)               = slosh_struct.f;
                amp_grid_all(indx_A,indx_T)             = slosh_struct.amp;
                A_end_grid_all(indx_A,indx_T)           = slosh_struct.A_end;
                T_grid_all(indx_A,indx_T)               = slosh_struct.T;
            end
        end
        
        % loading offset sloshing sensitivities and values
        for indx_A_offset = 1:N_A_offset
            load_string                                 = fullfile(location,sprintf('Sensitivity_offset_A_%.3f.mat',A_offset_list(indx_A_offset)));
            % disp(load_string)
            load(load_string,"slosh_struct")
            sensitivity_offset_all(indx_A_offset,:)        = slosh_struct.S_sloshy(1,:);
            f_offset(indx_A_offset)                        = slosh_struct.f;
        end
        
        % loading multiple linear ramp sensitivities and values
        for indx_A_ML = 1:N_A_ML
            for indx_T_ML = 1:N_T_ML
                load_string                                     = fullfile(location,sprintf('Sensitivity_CLR_A_end_%.3f_T_%.3f.mat',A_end_list_ML(indx_A_ML),T1_list_ML(indx_T_ML)+T2_list_ML(indx_T_ML)));
                % disp(load_string)
                load(load_string,"slosh_struct")
                sensitivity_f_ML_all(indx_A_ML,indx_T_ML,:)     = slosh_struct.S_sloshy(1,:);
                sensitivity_amp_ML_all(indx_A_ML,indx_T_ML,:)   = slosh_struct.S_sloshy(3,:);
                sensitivity_ph_ML_all(indx_A_ML,indx_T_ML,:)    = slosh_struct.S_sloshy(4,:);
                f_ML(indx_A_ML,indx_T_ML)                       = slosh_struct.f;
                amp_ML(indx_A_ML,indx_T_ML)                     = slosh_struct.amp;
                ph_ML(indx_A_ML,indx_T_ML)                      = slosh_struct.phase;
                A_end_grid_ML(indx_A_ML,indx_T_ML)              = slosh_struct.A_end;
                T_grid_ML(indx_A_ML,indx_T_ML)                  = slosh_struct.T;
            end
        end
        
        fprintf('loading data done!\n')
    else
        %% recalculate Sensitivities
        fprintf('-------------------------------------------------------------\nSensitivities are calculated new, this will take a while >3h!\nUse the heuristic method if you want a quick result!\n-------------------------------------------------------------\n')
        x0                  = [k1_est,k2_est,c_est,As_est,g1D_est];
        
        x_bounds            = 2.5;
        a6_coeffs           = 0;
        a4_coeffs           = @(x_is) x_is(2)/8/x_is(3)^2;
        a2_coeffs_pre_As    = @(x_is)[x_is(1)/2, -x_is(1)/2*x_is(4)];
        a2_coeffs_post_As   = @(x_is)[- x_is(2)/4, x_is(2)/4*x_is(4)];
        a2_fun              = @(x_is,A) polyval(heaviside(x_is(4)-A)*a2_coeffs_pre_As(x_is) + heaviside(A-x_is(4))*a2_coeffs_post_As(x_is),A);
        V_fun               = @(x_is,A,x) polyval([a4_coeffs(x_is),0,a2_fun(x_is,A),0,0],x).*heaviside(x+x_bounds).*heaviside(x_bounds-x)+1e5.*heaviside(-x-x_bounds)+1e5.*heaviside(-x_bounds+x);
        V_current           = @(A)V_fun(x0,A,x);
        par.g               = g1D_est;

        u0                  = calculate_groundstate_imag_time_GPE_ss(V_current(A_ramps_start),par,ones(size(x)),x,1e-5);
        %% Linear ramp sensitivities
        for indx_A = 1:N_A
            for  indx_T  = 1:N_T
                indx_total              = (indx_A-1)*N_T+indx_T;
                T                       = T_list(indx_T);
                A_t                     = linspace(A_0,A_end_list(indx_A),T/dt+1);
                tic
                sloshy                  = freq_and_amp_from_ramp(A_t,dt,V_current,x,par,u0,3);
                S_sloshy                = calculate_Sensitivity_Si_linear_ramp(A_0,A_end_list(indx_A),T,dt,V_fun,x0,x,par,u0);
                calc_time               = toc();
                slosh_struct.f          = sloshy(1);
                slosh_struct.amp        = sloshy(2);
                slosh_struct.d          = sloshy(3);
                slosh_struct.ph         = sloshy(4);
                slosh_struct.A_0        = A_0;
                slosh_struct.A_end      = A_end_list(indx_A);
                slosh_struct.T          = T_list(indx_T);
                slosh_struct.S_sloshy   = S_sloshy;
                sensitivity_f_all(indx_A,indx_T,:)      = S_sloshy(1,:);
                sensitivity_amp_all(indx_A,indx_T,:)    = S_sloshy(2,:);
                sensitivity_d_all(indx_A,indx_T,:)      = S_sloshy(3,:);
                sensitivity_ph_all(indx_A,indx_T,:)     = S_sloshy(4,:);
                f_grid_all(indx_A,indx_T)               = slosh_struct.f;
                amp_grid_all(indx_A,indx_T)             = slosh_struct.amp;
                A_end_grid_all(indx_A,indx_T)           = slosh_struct.A_end;
                T_grid_all(indx_A,indx_T)               = slosh_struct.T;
                save_string             = fullfile(location,sprintf('Sensitivity_LR_A_end_%.3f_T_%.3f.mat',A_end_list(indx_A),T_list(indx_T)));
                fprintf('(%i/%i)\n',indx_total,N_A*N_T)
                disp(strcat(save_string,sprintf(' in %.3f seconds',calc_time)))
                save(save_string,"slosh_struct")
    
            end
        end
        %% Consecutive Linear Ramps
        for indx_A_ML = 1:N_A_ML
            for  indx_T_ML  = 1:N_T_ML
                T1                      = T1_list_ML(indx_T_ML);
                T2                      = T2_list_ML(indx_T_ML);
                A_t                     = [linspace(A_0,A_1_list_ML(indx_A_ML),T1/dt+1),linspace(A_1_list_ML(indx_A_ML),A_end_list(indx_A_ML),T2/dt)];
                tic
                sloshy                  = freq_and_amp_from_ramp(A_t,dt,V_current,x,par,u0,3);
                S_sloshy                = calculate_Sensitivity_Si_multi_linear_ramp(A_0,A_1_list_ML(indx_A_ML),A_end_list(indx_A_ML),T1,T2,dt,V_fun,x0,x,par,u0);
                calc_time               = toc();
                slosh_struct.f          = sloshy(1);
                slosh_struct.amp        = sloshy(2);
                slosh_struct.d          = sloshy(3);
                slosh_struct.phase      = sloshy(4);
                slosh_struct.A_0        = A_0;
                slosh_struct.A_1        = A_1_list_ML(indx_A_ML);
                slosh_struct.A_end      = A_end_list(indx_A_ML);
                slosh_struct.T1         = T1_list_ML(indx_T_ML);
                slosh_struct.T2         = T2_list_ML(indx_T_ML);
                slosh_struct.T          = T1_list_ML(indx_T_ML)+T2_list_ML(indx_T_ML);
                slosh_struct.S_sloshy   = S_sloshy;
                sensitivity_f_ML_all(indx_A_ML,indx_T_ML,:)     = slosh_struct.S_sloshy(1,:);
                sensitivity_amp_ML_all(indx_A_ML,indx_T_ML,:)   = slosh_struct.S_sloshy(3,:);
                sensitivity_ph_ML_all(indx_A_ML,indx_T_ML,:)    = slosh_struct.S_sloshy(4,:);
                f_ML(indx_A_ML,indx_T_ML)                       = slosh_struct.f;
                amp_ML(indx_A_ML,indx_T_ML)                     = slosh_struct.amp;
                ph_ML(indx_A_ML,indx_T_ML)                      = slosh_struct.phase;
                A_end_grid_ML(indx_A_ML,indx_T_ML)              = slosh_struct.A_end;
                T_grid_ML(indx_A_ML,indx_T_ML)                  = slosh_struct.T;
                save_string             = fullfile(location,sprintf('Sensitivity_CLR_A_end_%.3f_T_%.3f.mat',A_end_list_ML(indx_A_ML),T1_list_ML(indx_T_ML)+T2_list_ML(indx_T_ML)));
                disp(strcat(save_string,sprintf(' in %.3f seconds',calc_time)))
                save(save_string,"slosh_struct")
    
            end
            fprintf('(%i/%i)\n',indx_A_ML,N_A_ML)
        end
        %% Offset Measurements
        for indx_A_offset = 1:N_A_offset
            A_i                     = A_offset_list(indx_A_offset);
            tic
            if indx_A_offset == 1
                u0_off                  = calculate_groundstate_imag_time_GPE_ss(V_current(A_i),par,u0,x,1e-5);
            else
                u0_off                  = calculate_groundstate_imag_time_GPE_ss(V_current(A_i),par,u0_off,x,1e-5);
            end
            [~,~,x0_freq]           = calculate_shaking_frequencies(A_i,V_current,x,par,1,u0_off,dt);
            S_sloshy                = calculate_Sensitivity_Si_offset(A_i,dt,V_fun,x0,x,par,u0_off);
            calc_time               = toc();
            slosh_struct.f          = x0_freq(1);
            slosh_struct.A          = A_offset_list(indx_A_offset);
            slosh_struct.S_sloshy   = S_sloshy;
            sensitivity_offset_all(indx_A_offset,:)        = slosh_struct.S_sloshy(1,:);
            f_offset(indx_A_offset)                        = slosh_struct.f;
            save_string             = fullfile(location,sprintf('Sensitivity_offset_A_%.3f.mat',A_offset_list(indx_A_offset)));
            disp(strcat(save_string,sprintf(' in %.3f seconds',calc_time)))
            save(save_string,"slosh_struct")
            fprintf('(%i/%i)\n',indx_A_offset,N_A_offset)
        end
        load(path_to_parameter,'data')
        data.is_done = true;
        save(path_to_parameter,'data')
    end

    %% run Genetic Algorithm
    run GA_internal
end



%% save results
N_linear_ramps              = numel(A_linear_ramps);
N_consecutive_linear_ramps  = numel(A_consecutive_linear_ramps);
N_offset                    = numel(A_offset);

% check for heuristic method
if use_heuristic_method
    calibration_experiment_struct.optimized_with_GA = false;
else
    calibration_experiment_struct.optimized_with_GA = true;
end

% linear ramps
calibration_experiment_struct.ramps_cell = {};
for indx = 1:N_linear_ramps
    t   = 0:dt:T_linear_ramps(indx);
    A_t = linspace(A_ramps_start,A_linear_ramps(indx),numel(t));
    ramps_struct.is_ramp   = 1;
    ramps_struct.A_0       = A_t(1);
    ramps_struct.A_end     = A_t(end);
    ramps_struct.T_ramp    = t(end);
    ramps_struct.t         = t;
    ramps_struct.A_t       = A_t;
    calibration_experiment_struct.ramps_cell{indx} = ramps_struct;
end

% consecutive linear ramps
for indx = 1:N_consecutive_linear_ramps
    t                       = 0:dt:6;
    t_linear                = 0:dt:1;
    t_after                 = t_linear(end)+dt:dt:t(end);
    A_t                     = [linspace(A_ramps_start,A_ramp_min,numel(t_linear)),linspace(A_ramp_min,A_consecutive_linear_ramps(indx),numel(t_after))];
    ramps_struct.is_ramp    = 1;
    ramps_struct.A_0        = A_t(1);
    ramps_struct.A_end      = A_t(end);
    ramps_struct.T_ramp     = t(end);
    ramps_struct.t          = t;
    ramps_struct.A_t        = A_t;
    calibration_experiment_struct.ramps_cell{indx+N_linear_ramps} = ramps_struct;
end

% offset measurements
calibration_experiment_struct.offset_cell   = {};
for indx = 1:N_offset
    offset_struct.is_ramp   = 0;
    offset_struct.A_offset  = A_offset(indx);
    calibration_experiment_struct.offset_cell{indx} = offset_struct;
end
calibration_experiment_struct.parameter_data    = data;


%% display resulting ramps
figure
for indx_ramp = 1:(numel(calibration_experiment_struct.ramps_cell))
    plot(calibration_experiment_struct.ramps_cell{indx_ramp}.t,calibration_experiment_struct.ramps_cell{indx_ramp}.A_t,'DisplayName',sprintf('%i: A = %.3f, T = %.3f',indx_ramp,calibration_experiment_struct.ramps_cell{indx_ramp}.A_end,calibration_experiment_struct.ramps_cell{indx_ramp}.t(end)))
    hold on
end
legend
grid on
legend('Location','northwest')


%% display results as grid map
[A_possible_linear_ramps_grid,T_possible_linear_ramps_grid] = meshgrid(A_possible_linear_ramps,T_possible_linear_ramps);

figure
subplot(5,1,[1,3])
plot(A_possible_linear_ramps_grid(:),T_possible_linear_ramps_grid(:),'.')
hold on
plot(A_linear_ramps,T_linear_ramps,'o','LineWidth',3,'MarkerSize',8)
xlim([A_offset_min-0.005,A_ramp_max+0.005])
grid on
title('Linear Ramps')
subplot(5,1,4)
plot(A_possible_CLR(:),ones(size(A_possible_CLR(:))),'.')
hold on
plot(A_consecutive_linear_ramps,ones(size(A_consecutive_linear_ramps)),'o','LineWidth',3,'MarkerSize',8)
xlim([A_offset_min-0.005,A_ramp_max+0.005])
grid on
title('Conceutive Linear Ramps')
subplot(5,1,5)
plot(A_possible_Offset(:),ones(size(A_possible_Offset(:))),'.')
hold on
plot(A_offset,ones(size(A_offset)),'o','LineWidth',3,'MarkerSize',8)
xlim([A_offset_min-0.005,A_ramp_max+0.005])
grid on
title('Offset')


% create selection struct
selection_struct.A_linear_ramps                 = A_linear_ramps;
selection_struct.T_linear_ramps                 = T_linear_ramps;
selection_struct.A_offset                       = A_offset;
selection_struct.A_consecutive_linear_ramps     = A_consecutive_linear_ramps;
selection_struct.parameter_data                 = data;