%% This is an internal function, be very carefull when changing anything!!!!

if ~use_estimates_from_selection
    % use original parameters fitted
    k1_est                  = -1.350570830057047e+03;   % slope of curvature before As
    k2_est                  = 2.392634373293948e+03;    % slope of curvature after As
    c_est                   = 2.239720917436003;        % factor before sqrt function for minimum position c*sqrt(A-As)
    As_est                  = 0.367412136099471;        % parameter at which the double well starts forming
    g1D_est                 = 6.838794525033462;        % nonlinear interaction constant for the transversal direction
else
    % use parameters used for experiment selection
    k1_est                  = calibration_experiment_struct.parameter_data.k1;  % slope of curvature before As
    k2_est                  = calibration_experiment_struct.parameter_data.k2;  % slope of curvature after As
    c_est                   = calibration_experiment_struct.parameter_data.c;   % factor before sqrt function for minimum position x_m = c*sqrt(A-As)
    As_est                  = calibration_experiment_struct.parameter_data.As;  % parameter at which the double well starts forming
    g1D_est                 = calibration_experiment_struct.parameter_data.g1D; % nonlinear interaction constant for the transversal direction
end

x0_unscaled     = [k1_est,k2_est,c_est,As_est,g1D_est]; % initial guess of parameters

par_norm        = [1e3,1e3,2,0.4,6]; % norm values for the parameters
x0              = x0_unscaled./par_norm; % normalized initial guesses of parameters

% lower and upper bound for the parameters
lb              = [-3000,500,1.5,0.3,5]./par_norm; 
ub              = [-500,4000,3,0.5,8]./par_norm;

offset_shift    = 3*ones(size(f_values_offset_exp)); % values of initial shift of ground state in offset measurements in discretization steps


%% Defining potential
x_bounds            = 2.5; % max range of potential
a6_coeffs           = @(x_is)0;   % coefficients of a6(A)*x6 part
a4_coeffs           = @(x_is)x_is(2)/8/x_is(3)^2; % coefficients of a4(A)*x4 part
a2_coeffs_pre_As    = @(x_is)[x_is(1)/2, -x_is(1)/2*x_is(4)];  % coefficients of a2(A)*x2 part before As
a2_coeffs_post_As   = @(x_is)[- x_is(2)/4, x_is(2)/4*x_is(4)]; % coefficients of a2(A)*x2 part after As
a2_fun              = @(x_is,A) polyval(heaviside(x_is(4)-A)*a2_coeffs_pre_As(x_is) + heaviside(A-x_is(4))*a2_coeffs_post_As(x_is),A);  % combined function of a2(A)*x^2
V_fun               = @(x_is,A,x) polyval([a4_coeffs(x_is),0,a2_fun(x_is,A),0,0],x).*heaviside(x+x_bounds).*heaviside(x_bounds-x)+1e5.*heaviside(-x-x_bounds)+1e5.*heaviside(-x_bounds+x); % combined V as function x and parameters
V_current           = @(A)V_fun(x0.*par_norm,A,x); % initial Potential as function of A
V_fun_scaled        = @(x_scaled,A,x) V_fun(x_scaled.*par_norm,A,x); % Potential as function of scaled parameters, A and x
par.g               = g1D_est; % initial g1D/g_perp


% extracting ramps and offset values
exp_ramps   = calibration_experiment_struct.ramps_cell;
exp_offset  = calibration_experiment_struct.offset_cell;
N_ramps     = size(exp_ramps,2);
N_offset    = size(exp_offset,2);

% check size of measurements
if N_ramps ~= size(f_values_ramps_exp,2)
    disp('Number of ramp measurements is not the same as number of ramps from provided calibration selection!!!!')
    return
end
if N_offset ~= size(f_values_offset_exp,2)
    disp('Number of offset measurements is not the same as number of offset parameters from provided calibration selection!!!!')
    return
end

% rescale ramps to correct time discretization
if exp_ramps{1}.t(2) - exp_ramps{1}.t(1) ~= dt
    for indx_ramp = 1:N_ramps
        t_new   = 0:dt:exp_ramps{indx_ramp}.T_ramp;
        A_t_new = interp1(exp_ramps{indx_ramp}.t,exp_ramps{indx_ramp}.A_t,t_new,"linear");
        exp_ramps{indx_ramp}.t      = t_new;
        exp_ramps{indx_ramp}.A_t    = A_t_new;
    end
end

%% initialization
% get initial parameter set
par_init    = par;
par_init.g  = x0(end)*par_norm(end);
V_inital    = @(A) V_fun_scaled(x0,A,x);

% get parameter data for initial point of ramps
parameter_data  = calibration_experiment_struct.parameter_data;
A_ramps_start   = parameter_data.A_ramps_start;
A_plot          = linspace(A_ramps_start,parameter_data.A_ramp_max,1e3); % range for plotting

fprintf('calculating initial states...\n')
% initial state for ramps
u0_init         = calculate_groundstate_imag_time_GPE_ss(V_inital(A_ramps_start),par_init,ones(size(x)),x);

% initial states for offset measurements
u0_offset_init          = zeros(N_x,N_offset);
for indx_offset = 1:N_offset
    if indx_offset == 1
        u0_offset_init(:,indx_offset)   = calculate_groundstate_imag_time_GPE_ss(V_inital(exp_offset{indx_offset}.A_offset),par_init,u0_init,x);
    else
        u0_offset_init(:,indx_offset)   = calculate_groundstate_imag_time_GPE_ss(V_inital(exp_offset{indx_offset}.A_offset),par_init,u0_offset_init(:,indx_offset-1),x);
    end
end
fprintf('done!\n')
fprintf('calculating initial values...\n')
% calculate initial values
f_ramps_init    = zeros(size(f_values_ramps_exp));
amp_ramps_init  = zeros(size(f_values_ramps_exp));
d_ramps_init    = zeros(size(f_values_ramps_exp));
ph_ramps_init   = zeros(size(f_values_ramps_exp));

f_offset_init   = zeros(size(f_values_offset_exp));

for indx_ramp = 1:N_ramps
    sloshy                      = freq_and_amp_from_ramp(exp_ramps{indx_ramp}.A_t,dt,V_inital,x,par_init,u0_init,3);
    f_ramps_init(1,indx_ramp)   = sloshy(1);
    amp_ramps_init(indx_ramp)   = sloshy(2);
    d_ramps_init(indx_ramp)     = sloshy(3);
    ph_ramps_init(indx_ramp)    = sloshy(4);
    if size(f_values_ramps_exp,1)>1
        if f_values_ramps_exp(1,indx_ramp) ~= f_values_ramps_exp(2,indx_ramp)
            f_ramps_init(:,indx_ramp)   = sort(sloshy([1,11]));
        else
            f_ramps_init(:,indx_ramp)   = sloshy(1);
        end
    else
        f_ramps_init(:,indx_ramp)   = sloshy(1);
    end
end

for indx_offset = 1:N_offset
    [~,~,x0_freq,x0_freq_both]  = calculate_shaking_frequencies(exp_offset{indx_offset}.A_offset,V_inital,x,par_init,1,u0_offset_init(:,indx_offset),dt,offset_shift(indx_offset));

    if numel(f_offset_init) == N_offset
        f_offset_init(indx_offset)  = x0_freq;
    else
        if f_values_offset_exp(1,indx_offset) == f_values_offset_exp(2,indx_offset)        
            f_offset_init(:,indx_offset) = [x0_freq;x0_freq];
        else
            f_offset_init(:,indx_offset) = x0_freq_both;
        end
    end
end
fprintf('done!\n')

% plotting comparisson of experiment and initial simulation
figure
subplot(2,2,1)
plot(1:N_ramps,f_values_ramps_exp,"Color",[0 0.4470 0.7410],'LineStyle','none','LineWidth',3,'Marker','x','MarkerSize',10,'DisplayName','Experiment')
hold on
plot(1:N_ramps,f_ramps_init,"Color",[0.8500 0.3250 0.0980],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
grid on
ylabel('Frequencies in kHz')
title('Ramps')
subplot(2,2,2)
plot(1:N_offset,f_values_offset_exp,"Color",[0 0.4470 0.7410],'LineStyle','none','LineWidth',3,'Marker','x','MarkerSize',10,'DisplayName','Experiment')
hold on
plot(1:N_offset,f_offset_init,"Color",[0.8500 0.3250 0.0980],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
grid on
legend
title('Offset')
subplot(2,2,3)
plot(1:N_ramps,abs(f_values_ramps_exp-f_ramps_init),"Color","black",'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
grid on
ylabel('Error of frequencies in kHz')
subplot(2,2,4)
plot(1:N_offset,abs(f_values_offset_exp-f_offset_init),"Color","black",'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
grid on
drawnow

%% define optimization problem
% condition parameters of optimization problem
A_bcond             = [];
b_bcond             = [];
Aeq                 = [];
Beq                 = [];
nonlinconst         = [];

% cost function
cost_fun    = @(params) cost_ramps_and_offset_only_f(params,V_fun_scaled,x,par,dt,exp_ramps,exp_offset,u0_init,u0_offset_init,f_values_ramps_exp,f_values_offset_exp,par_norm,offset_shift);

% options for cost function
opti_opts   = optimoptions('fmincon','Display','iter-detailed','SpecifyConstraintGradient',false,'EnableFeasibilityMode',false,'Algorithm','interior-point','StepTolerance',1e-8,'FiniteDifferenceStepSize',1e-3,'UseParallel',true);

% test cost function evaluation and time
fprintf('testing initial cost function evaluation\n')
tic
cost0   = cost_fun(x0);
T0      = toc();
fprintf('initial cost: %.5f in %.5f s\n',cost0,T0)

% set further options here for easy availabillity
opti_opts.MaxFunctionEvaluations    = 1e4;
opti_opts.MaxIterations             = 1e4;
fprintf('starting optimization\n')


%% Run optimization
tic
[params_final,fval,exitflag,output,lambda,grad,hessian]   = fmincon(cost_fun,x0,A_bcond,b_bcond,Aeq,Beq,lb,ub,nonlinconst,opti_opts);
T_opt = toc();
fprintf('optimization done in %.3f seconds\n',T_opt)

%% extract results from optimization
% check if lower than final fval was found and use that
if output.bestfeasible.fval < fval
    params_opt_scaled = output.bestfeasible.x;
    fval       = output.bestfeasible.fval;
else
    params_opt_scaled = params_final;
end

% set optimal parameters with calibrated g1D
params_opt      = params_opt_scaled.*par_norm;
par_opt    = par;
par_opt.g = params_opt_scaled(end)*par_norm(end);

V_opt           = @(A)V_fun(params_opt,A,x);

%% calculate results of calibration for comparisson
% calculate optimal initial state
u0_opt          = calculate_groundstate_imag_time_GPE_ss(V_opt(A_ramps_start),par_opt,u0_init,x);
% initial  states for offset measurements optimized
u0_offset_opt          = zeros(N_x,N_offset);
for indx_offset = 1:N_offset
    if indx_offset == 1
        u0_offset_opt(:,indx_offset)   = calculate_groundstate_imag_time_GPE_ss(V_opt(exp_offset{indx_offset}.A_offset),par_opt,u0_opt,x);
    else
        u0_offset_opt(:,indx_offset)   = calculate_groundstate_imag_time_GPE_ss(V_opt(exp_offset{indx_offset}.A_offset),par_opt,u0_offset_opt(:,indx_offset-1),x);
    end
end

% calculate initial values
f_ramps_opt    = zeros(size(f_values_ramps_exp));
amp_ramps_opt  = zeros(size(f_values_ramps_exp));
d_ramps_opt    = zeros(size(f_values_ramps_exp));
ph_ramps_opt   = zeros(size(f_values_ramps_exp));

f_offset_opt   = zeros(size(f_values_offset_exp));

for indx_ramp = 1:N_ramps
    sloshy                      = freq_and_amp_from_ramp(exp_ramps{indx_ramp}.A_t,dt,V_opt,x,par_opt,u0_init,3);
    f_ramps_opt(indx_ramp)     = sloshy(1);
    amp_ramps_opt(indx_ramp)   = sloshy(2);
    d_ramps_opt(indx_ramp)     = sloshy(3);
    ph_ramps_opt(indx_ramp)    = sloshy(4);
    if f_values_ramps_exp(1,indx_ramp) ~= f_values_ramps_exp(2,indx_ramp)
        f_ramps_opt(:,indx_ramp)    = sort(sloshy([1,11]));
    else
        f_ramps_opt(:,indx_ramp)    = sloshy(1);
    end
end

for indx_offset = 1:N_offset
    [~,~,x0_freq,x0_freq_both]  = calculate_shaking_frequencies(exp_offset{indx_offset}.A_offset,V_opt,x,par_opt,1,u0_offset_opt(:,indx_offset),dt,offset_shift(indx_offset));

    if numel(f_offset_opt) == N_offset
        f_offset_opt(indx_offset)  = x0_freq;
    else
        if f_values_offset_exp(1,indx_offset) == f_values_offset_exp(2,indx_offset)        
            f_offset_opt(:,indx_offset) = [x0_freq;x0_freq];
        else
            f_offset_opt(:,indx_offset) = x0_freq_both;
        end
    end
end

%% extracting model parameters and polynomail coefficients from optimized parameters vector
k1_opt                  = params_opt(1);
k2_opt                  = params_opt(2);
c_opt                   = params_opt(3);
As_opt                  = params_opt(4);
g1D_opt                 = params_opt(5);
coeffs_a6_opt           = 0;
coeffs_a4_opt           = k2_opt/8/c_opt^2;
coeffs_a2_pre_As_opt    = [k1_opt/2, -k1_opt/2*As_opt];
coeffs_a2_post_As_opt   = [- k2_opt/4, k2_opt/4*As_opt];
coeffs_a0_opt           = 0;


%% plotting comparisson of experiment and initial and optimized simulation
figure
subplot(2,2,1)
plot(1:N_ramps,f_values_ramps_exp,"Color",[0 0.4470 0.7410],'LineStyle','none','LineWidth',3,'Marker','x','MarkerSize',10,'DisplayName','Experiment')
hold on
plot(1:N_ramps,f_ramps_init,"Color",[0.8500 0.3250 0.0980],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
plot(1:N_ramps,f_ramps_opt,"Color",[0.9290 0.6940 0.1250],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Calibrated')
grid on
legend
ylabel('Frequencies in kHz')
title('Ramps')
subplot(2,2,2)
plot(1:N_offset,f_values_offset_exp,"Color",[0 0.4470 0.7410],'LineStyle','none','LineWidth',3,'Marker','x','MarkerSize',10,'DisplayName','Experiment')
hold on
plot(1:N_offset,f_offset_init,"Color",[0.8500 0.3250 0.0980],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
plot(1:N_offset,f_offset_opt,"Color",[0.9290 0.6940 0.1250],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Calibrated')
grid on
legend
title('Offset')
subplot(2,2,3)
plot(1:N_ramps,abs(f_values_ramps_exp-f_ramps_init).^2,"Color",[0.8500 0.3250 0.0980],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
hold on
plot(1:N_ramps,abs(f_values_ramps_exp-f_ramps_opt).^2,"Color",[0.9290 0.6940 0.1250],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Calibrated')
grid on
legend
subplot(2,2,4)
plot(1:N_offset,abs(f_values_offset_exp-f_offset_init).^2,"Color","black",'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Initial')
hold on
plot(1:N_offset,abs(f_values_offset_exp-f_offset_opt).^2,"Color",[0.9290 0.6940 0.1250],'LineStyle','none','LineWidth',3,'Marker','o','MarkerSize',10,'DisplayName','Calibrated')
grid on
legend
