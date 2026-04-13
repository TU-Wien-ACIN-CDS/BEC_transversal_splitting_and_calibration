%% This is an internal function, be very carefull when changing anything!!!!

%% parameters of optimization
maxIterBFGS             = 350;  % max iterations per scenario
n_old_max               = 8;   % max iteration used for BFGS approximation
max_NE_steps            = 15;   % max amount of not effective iterations allowed before scenario is abandoned
choose_cost             = 1;    % if 1: use energy cost, other options are not tested
par_bfgs.use_memLimit   = 1;    % if 1, use memory limited BFGS
par.gamma               = 2e-3;  % regularization factor for lambda cost
par_bfgs.A_min          = 0;    % minimum value of A_t allowed
par_bfgs.A_max          = 1;    % maximum value of A_t allowed
par_bfgs.use_H1_cost    = 1;    % if 1, use H1 cost (see "Computational techniques for a quantum control problem with H1-cost"), if 0, use L2 cost, if 2, use combined H1 and L2 cost, if 3 use smoothing of L2 cost
par_bfgs.C_H1_cost      = 1e-5; % scaling factor of L2 component of combined cost, only used if par_bfgs.use_H1_cost == 2, standard value = 1e-5

%% Defining potential
parameters_calibrated   = [calib_struct.k1,calib_struct.k2,calib_struct.c,calib_struct.As];
par.g                   = calib_struct.g1D*0;

g1D                 = par.g;
m                   = par.m;


x_bound             = 2.5; % max range of potential
a6_coeffs           = 0;  % coefficients of a6(A)*x6 part
a4_coeffs           = @(x_is)x_is(2)/8/x_is(3)^2; % coefficients of a4(A)*x4 part
a2_coeffs_pre_As    = @(x_is)[x_is(1)/2, -x_is(1)/2*x_is(4)];  % coefficients of a2(A)*x2 part before As
a2_coeffs_post_As   = @(x_is)[- x_is(2)/4, x_is(2)/4*x_is(4)]; % coefficients of a2(A)*x2 part after As
a2_fun              = @(x_is,A) polyval(heaviside(x_is(4)-A)*a2_coeffs_pre_As(x_is) + heaviside(A-x_is(4))*a2_coeffs_post_As(x_is),A);  % combined function of a2(A)*x^2
V_fun               = @(x_is,A,x) polyval([a4_coeffs(x_is),0,a2_fun(x_is,A),0,0],x).*heaviside(x+x_bound).*heaviside(x_bound-x)+1e5.*heaviside(-x-x_bound)+1e5.*heaviside(-x_bound+x); % combined V as function x and parameters
V_A                 = @(A) V_fun(parameters_calibrated,A,x); % Potential as function of A

da6_poly_calib     = 0; % coefficients of derivative of a6(A)
da4_poly_calib     = 0; % coefficients of derivative of a4(A)
da2_poly_pre_As    = polyder(a2_coeffs_pre_As(parameters_calibrated));  % coefficients of derivative of a2(A) before As
da2_poly_post_As   = polyder(a2_coeffs_post_As(parameters_calibrated));  % coefficients of derivative of a2(A) after As
delta_A            = 0.07;  % range of A value for smothing of a2(A) derivative around As
da_inter_fun       = @(A) heaviside(-calib_struct.As+delta_A+A).*heaviside(-A+delta_A+calib_struct.As).*((polyval(da2_poly_post_As,calib_struct.As+delta_A)-polyval(da2_poly_pre_As,calib_struct.As-delta_A)).*(1./(1+exp(-6/delta_A.*(A-calib_struct.As))))+polyval(da2_poly_pre_As,calib_struct.As-delta_A)); % function used for smoothing
da2_fun            = @(A) heaviside(calib_struct.As-delta_A-A).*polyval(da2_poly_pre_As,A)+heaviside(A-delta_A-calib_struct.As).*polyval(da2_poly_post_As,A)+da_inter_fun(A); % function of derivative of a2
da0_poly_calib     = 0; % derivative of a0(A)
dV_A               = @(A,n) (polyval([polyval(da6_poly_calib,A),0,polyval(da4_poly_calib,A),0,da2_fun(A),0,polyval(da0_poly_calib,A)],x).*heaviside(x+x_bound).*heaviside(x_bound-x)); % derivative of Potnetial
    

ham                 = @(A)-lap4/(2*m)+spdiag(V_A(A)); % hamiltonian component
energy_fun          = @(u,A) real(trapz(x,conj(u).*(ham(A)+par.g/2*(spdiag(abs(u).^2)))*u)); % function to calculate energy
lambda_cost         = @(t,lambda_t) par.gamma * mean(trapz(t(1:end-1),diff(lambda_t,1).^2,2),1); % function to calculate cost of regularization component


%% design scenario
u0_case             = zeros(N_x,N_ramps);
u1_case             = zeros(N_x,N_ramps);
t_case_ramp         = cell(1,N_ramps);
N_t_case            = zeros(1,N_ramps);
A_t_init_case       = cell(1,N_ramps);
A_t_init_used       = cell(1,N_ramps);
A_t_is_cases        = cell(1,N_ramps);
t_case              = cell(1,N_ramps);

fprintf('Calculating initial and final ground states...\n')
switch switch_scenario
    case 1 % linear ramps for set A_0, A_end and T values
        for indx_ramp = 1:N_ramps
            t_case_ramp{indx_ramp}            = 0:dt:T_list(indx_ramp);
            N_t_case(indx_ramp)          = numel(t_case_ramp{indx_ramp});
            
            A_t_init_case{indx_ramp}     = linspace(A_0_list(indx_ramp),A_end_list(indx_ramp),N_t_case(indx_ramp));

            if indx_ramp > 1
                if A_0_list(indx_ramp) ~= A_0_list(indx_ramp-1)
                    u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_0_list(indx_ramp)),par,ones(size(x)),x);
                else
                    u0_case(:,indx_ramp)        = u0_case(:,indx_ramp-1);
                end
                if A_end_list(indx_ramp) ~= A_end_list(indx_ramp-1)
                    u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,ones(size(x)),x);
                else
                    u1_case(:,indx_ramp)        = u1_case(:,indx_ramp-1);
                end
            else
                u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_0_list(indx_ramp)),par,ones(size(x)),x);
                u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,u0_case(:,indx_ramp),x);
            end
        end
        
    case 2 % linear ramps for A_0, A_end and T values set by previously run scenario
        T_list              = zeros(1,N_ramps);
        A_0_list            = zeros(1,N_ramps);
        A_end_list          = zeros(1,N_ramps);
        
        for indx_ramp = 1:N_ramps
            A_0_list(indx_ramp)     = optimized_trajectories_struct.ramps_cell{indx_ramp}.A_0;
            A_end_list(indx_ramp)   = optimized_trajectories_struct.ramps_cell{indx_ramp}.A_end;
            T_list(indx_ramp)       = optimized_trajectories_struct.ramps_cell{indx_ramp}.T;

            t_case_ramp{indx_ramp}            = 0:dt:T_list(indx_ramp);
            N_t_case(indx_ramp)          = numel(t_case_ramp{indx_ramp});
            
            A_t_init_case{indx_ramp}     = linspace(A_0_list(indx_ramp),A_end_list(indx_ramp),N_t_case(indx_ramp));

            if indx_ramp > 1
                if A_0_list(indx_ramp) ~= A_0_list(indx_ramp-1)
                    u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_0_list(indx_ramp)),par,ones(size(x)),x);
                else
                    u0_case(:,indx_ramp)        = u0_case(:,indx_ramp-1);
                end
                if A_end_list(indx_ramp) ~= A_end_list(indx_ramp-1)
                    u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,ones(size(x)),x);
                else
                    u1_case(:,indx_ramp)        = u1_case(:,indx_ramp-1);
                end
            else
                u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_0_list(indx_ramp)),par,ones(size(x)),x);
                u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,ones(size(x)),x);
            end
        end

    case 3 % optimized ramps from previously run scenario
        T_list              = zeros(1,N_ramps);
        A_0_list            = zeros(1,N_ramps);
        A_end_list          = zeros(1,N_ramps);
        
        for indx_ramp = 1:N_ramps
            A_0_list(indx_ramp)     = optimized_trajectories_struct.ramps_cell{indx_ramp}.A_0;
            A_end_list(indx_ramp)   = optimized_trajectories_struct.ramps_cell{indx_ramp}.A_end;
            T_list(indx_ramp)       = optimized_trajectories_struct.ramps_cell{indx_ramp}.T;

            t_case_ramp{indx_ramp}            = 0:dt:T_list(indx_ramp);
            N_t_case(indx_ramp)          = numel(t_case_ramp{indx_ramp});
            
            A_t_init_case{indx_ramp}     = interp1(optimized_trajectories_struct.ramps_cell{indx_ramp}.t,optimized_trajectories_struct.ramps_cell{indx_ramp}.A_t,t_case_ramp{indx_ramp}, 'nearest','extrap');

            if indx_ramp > 1
                if A_0_list(indx_ramp) ~= A_0_list(indx_ramp-1)
                    u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_0_list(indx_ramp)),par,ones(size(x)),x);
                else
                    u0_case(:,indx_ramp)        = u0_case(:,indx_ramp-1);
                end
                if A_end_list(indx_ramp) ~= A_end_list(indx_ramp-1)
                    u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,ones(size(x)),x);
                else
                    u1_case(:,indx_ramp)        = u1_case(:,indx_ramp-1);
                end
            else
                u0_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_0_list(indx_ramp)),par,ones(size(x)),x);
                u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,ones(size(x)),x);
            end
        end

    case 4 % linear ramps starting from a state and A_0 from a simulation
        for indx_ramp = 1:N_ramps

            A_0_list(indx_ramp)         = initial_state_struct.A_initial;

            t_case_ramp{indx_ramp}           = 0:dt:T_list(indx_ramp);
            N_t_case(indx_ramp)         = numel(t_case_ramp{indx_ramp});
            
            A_t_init_case{indx_ramp}    = linspace(A_0_list(indx_ramp),A_end_list(indx_ramp),N_t_case(indx_ramp));

            u0_case(:,indx_ramp)        = initial_state_struct.Psi_end;

            if indx_ramp > 1
                if A_end_list(indx_ramp) ~= A_end_list(indx_ramp-1)
                    u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,u1_case(:,indx_ramp-1),x);
                else
                    u1_case(:,indx_ramp)        = u1_case(:,indx_ramp-1);
                end
            else
                u1_case(:,indx_ramp)        = calculate_groundstate_imag_time_GPE_ss(V_A(A_end_list(indx_ramp)),par,ones(size(x)),x);
            end
        end

end
fprintf('Done\n')

%% starting optimization
for indx_ramp = 1:N_ramps
    
    figure(12)
    set(gcf,'Name','Optimized trajectory')
    title(sprintf('IOA %i out of %i: trajectory',indx_ramp,N_ramps))
    % hold off

    figure(11)
    set(gcf,'Name','Current Iterations')
    
    % get parameters for current optimization
    N_t         = N_t_case(indx_ramp);
    T           = T_list(indx_ramp);
    t           = t_case_ramp{indx_ramp};
    A0          = A_0_list(indx_ramp);
    A_end       = A_end_list(indx_ramp);
    A_t_init    = A_t_init_case{indx_ramp};
    u0          = u0_case(:,indx_ramp);
    u1          = u1_case(:,indx_ramp);

    % select if previous ramp is reused
    if and(switch_reuse==1,indx_ramp > 1)
        N_t_old     = numel(A_t_is);
        A_t_is_temp = smooth(interp1(linspace(0,1,N_t_old),A_t_is,linspace(0,1,N_t),"spline"),'loess')';
        A_t_is      = (A_t_is_temp - A_t_is_temp(1))*(A_end - A0)/(A_t_is_temp(end) - A_t_is_temp(1)) + A0;
    else
        A_t_is      = A_t_init;
    end
    A_t_init_used{indx_ramp}    = A_t_is;
    
    % function to calculate second derivative in time
    par_bfgs.calc_second_deriv_t = @(a)calculate_second_deriv(a,finit_difference_matrix_norm(3,:),n_order,t(2)-t(1));

    % storage for results
    cost_vals               = zeros(1,maxIterBFGS);
    lambda_cost_val         = zeros(1,maxIterBFGS);
    A_t_vals                = zeros(N_t,maxIterBFGS);
    a_opt_vals              = zeros(1,maxIterBFGS);
    grad_vals               = A_t_vals;
    search_direction_vals   = grad_vals;
    iter_times              = cost_vals;
    iter_n_iters_opt        = cost_vals;
    Pop_distribution_t      = zeros(2,maxIterBFGS);

    state_cost  = @(u) 1/2*(1-abs(trapz(x,conj(u1).*u)).^2); % function to calculate state cost
    mue_0       = par.g * max(abs(u0)).^2; % chemical potential for initial potential
    mue_1       = par.g * max(abs(u1)).^2; % chemical potential for final potential
    
    switch choose_cost
        case 0 % state cost
            end_val_fun         = @(u,A_end) 1i*trapz(x,u1.*u).*u1;
            cost_fun            = @(psi,lambda_t,t) state_cost(psi(:,end)) + lambda_cost(t,lambda_t);
            cost_fun_barebones  = @(lambda_t) cost_val_multi(cost_fun,N_t,x,u0,dt,V_A,lambda_t,par);
            % val_opt_thresh      = 1e-5;
        case 1 % energy cost
            end_val_fun         = @(u,lambda_end) -2i*(ham(lambda_end)+spdiag(par.g*abs(u).^2))*u;
            cost_fun            = @(psi,lambda_t,t) energy_fun(psi(:,end),lambda_t(:,end)) - energy_fun(u1,lambda_t(:,end)) + lambda_cost(t,lambda_t);
            cost_fun_barebones  = @(lambda_t) cost_val_multi(cost_fun,N_t,x,u0,dt,V_A,lambda_t,par);
            % val_opt_thresh      = 1e-4;
        case 2 % alpha * state cost + beta  * energy cost
            alpha               = 0.1;
            beta                = 1;
            end_val_fun         = @(u,A_end) alpha * (1i*trapz(x,conj(u1).*u).*u1) + beta * (-2i*(ham(A_end)+spdiag(par.g*abs(u).^2))*u) * 2 * (energy_fun(u,A_end) - energy_fun(u1,A_end)); 
            cost_fun            = @(psi,A_t,t) alpha * state_cost(psi(:,end)) + beta * (energy_fun(psi(:,end),A_t(:,end)) - energy_fun(u1,A_t(:,end)))^2 + lambda_cost(t,A_t);
            cost_fun_barebones   = @(A_t) cost_val_multi(cost_fun,N_t,x,u0,dt,V_A,A_t,par);
        case 3 % diff in energy ^2 cost
            end_val_fun         = @(u,A_end) (-2i*(ham(A_end)+spdiag(par.g*abs(u).^2))*u) * 2 * (energy_fun(u,A_end) - energy_fun(u1,A_end)); 
            % e_u1                = ;
            cost_fun            = @(psi,A_t,t) (energy_fun(psi(:,end),A_t(:,end)) - energy_fun(u1,A_t(:,end)))^2 + lambda_cost(t,A_t);
            cost_fun_barebones  = @(A_t) cost_val_multi(cost_fun,N_t,x,u0,dt,V_A,A_t,par);
            cost_fun            = @(psi,A_t,t) (energy_fun(psi(:,end),A_t(:,end)) - energy_fun(u1,A_t(:,end)))^2 + lambda_cost(t,A_t);
            cost_fun_barebones  = @(A_t) cost_val_multi(cost_fun,N_t,x,u0,dt,V_A,A_t,par);
    end
    
    % collect constant paramters for BFGS
    par_bfgs.n_old_max                  = n_old_max;
    lap_t                               = sparse(imfilter(eye(N_t),finit_difference_matrix_t(3,:)));
    par_bfgs.gradient_second_deriv      = @(psi_store,p_store,dV_lambda,lambda,second_deriv_old) calculate_J_lambda_pp(psi_store,p_store,dV_lambda,lambda,second_deriv_old,x,par);
    par_bfgs.lap_t                      = lap_t;
    par_bfgs.use_numerical_derivative   = 0;
    par_bfgs.dt                         = dt;
    
    % initial parameters for BFGS
    n_old           = 0;
    ss_old          = zeros(N_t,n_old_max,1);
    ys_old          = zeros(N_t,n_old_max,1);
    gradients_old   = zeros(N_t,1);
    n_lambdas       = 1;
    NE_steps        = 0;
    x_bound         = 2.5;
    x_selector      = abs(x) < x_bound;

    % plot of 
    figure(11)
    subplot(4,1,1)
    hold off
    plot(t,A_t_is)
    grid on
    hold off
    set(gca,'Xticklabel',[]);
    title(sprintf('Ramp %i out of %i: search iteration %i',indx_ramp,N_ramps,BFGS_iter))
    xlim([t(1),t(end)])
    xlabel('t in ms')
    ylabel('A')
    subplot(4,1,2)
    hold off
    plot(t,zeros(size(t)))
    xlim([t(1),t(end)])
    xlabel('t in ms')
    ylabel('delta A')
    title('search direction')
    subplot(4,1,3)
    hold off
    plot(t,zeros(size(t)))
    grid on
    hold off
    title('interwell distance')
    xlim([t(1),t(end)])
    xlabel('t in ms')
    ylabel('d in µm')
    subplot(4,1,4)
    hold off
    imagesc(t,x,0)
    title('Density')
    xlabel('t in ms')
    ylabel('x in µm')
    drawnow
    

    Hk = 1*eye(N_t); % matrix for non memory limited BFGS
    if ~par_bfgs.use_memLimit
        n_old = Hk;
    end

    % old derivative
    second_derivs_old               = smooth(calculate_second_deriv(A_t_is,finit_difference_matrix_norm(3,:),n_order,dt),0.1,'rloess')'; 

    for BFGS_iter = 1:maxIterBFGS
        fprintf('starting step %i of optimization %i / %i\n',BFGS_iter,indx_ramp,N_ramps)
        tic
        [psi_store]                 = calc_prop_GPE(N_t,x,u0,dt,V_A,A_t_is,par); % calculate state for current trajectory
        current_energy_diff         = energy_fun(psi_store(:,end),A_t_is(end))-energy_fun(u1,A_t_is(end)); % calculate current energy cost
        [interwell_dist,sigma_traj] = calculate_difference_from_fit(x,psi_store); % calculate interwell distance and sigma for current trajectory
        cost_0 = cost_fun(psi_store,A_t_is,t);
        if BFGS_iter == 1           % if initial
            interwell_dist_0    = interwell_dist;
            A_t_0               = A_t_is;
            initial_cost        = cost_0;
        end
        fprintf('with energy cost %.5e\n',current_energy_diff)

        if  cost_0 < val_opt_thresh  && BFGS_iter>1 % if end condition is met
            energy_end(indx_ramp)           = energy_fun(psi_store(:,end),A_t_is(end));
            cost_vals(BFGS_iter:end)        = cost_fun_barebones(A_t_is);
            lambda_cost_val(BFGS_iter:end)  = lambda_cost(t,A_t_is);
            fprintf('threshold reached\n')
            break;
        end
        
        % propagation of p
        % p_prop_function             = @(psi_store,lambda_t)calc_prop_GPE_p_v2(psi_store, lambda_t, x, lap4, V_A, m, g1D, mue_0, t, end_val_fun, x_selector);
        p_prop_function             = @(psi_store,lambda_t)calc_prop_GPE_p_v3(psi_store, lambda_t, x, lap_scaled, ilap_scaled, V_A, g1D, mue_0, t, end_val_fun);

        % plotting current state and trajectory
        figure(11)
        subplot(4,1,1)
        plot(t,A_t_is)
        grid on
        hold on
        plot(t,A_t_0,'k--')
        hold off
        set(gca,'Xticklabel',[]);
        title(sprintf('Ramp %i out of %i: search iteration %i',indx_ramp,N_ramps,BFGS_iter))
        xlim([t(1),t(end)])
        xlabel('t in ms')
        ylabel('A')
        subplot(4,1,2)
        xlim([t(1),t(end)])
        xlabel('t in ms')
        ylabel('delta A in µm')
        title('search direction')
        subplot(4,1,3)
        plot(t,interwell_dist)
        grid on
        hold on
        plot(t,interwell_dist_0,'k--')
        hold off
        title('interwell distance')
        ylim([0,min([max(interwell_dist)+0.1,5])])
        xlim([t(1),t(end)])
        xlabel('t in ms')
        ylabel('d in µm')
        subplot(4,1,4)
        hold off
        imagesc(t,x,abs(psi_store).^2)
        title('Density')
        xlabel('t in ms')
        ylabel('x in µm')
        drawnow

        % make BFGS step
        [A_t_is,ys_old,ss_old,n_old,gradients_old,search_directions,second_derivs_old,val_opt,a_opt,n_iter_opt,breaker,Hk] = make_BFGS_step(x,p_prop_function,psi_store,V_A,dV_A,A_t_is,ys_old,ss_old,gradients_old,second_derivs_old,n_old,n_lambdas,cost_fun_barebones,par_bfgs,Hk,cost_0);
        
        % store results of BFGS step
        lambda_cost_val(BFGS_iter)          = lambda_cost(t,A_t_is);
        cost_vals(BFGS_iter)                = val_opt;
        A_t_vals(:,BFGS_iter)               = A_t_is(:);
        grad_vals(:,BFGS_iter)              = gradients_old(:,end,1);   
        search_direction_vals(:,BFGS_iter)  = search_directions(:,end,1);                        
        iter_time                           = toc();
        iter_times(BFGS_iter)               = iter_time;
        iter_n_iters_opt(BFGS_iter)         = n_iter_opt;
        N1_t                                = trapz(x(1:ceil(N_x/2)),abs(psi_store(1:ceil(N_x/2))).^2);
        N2_t                                = trapz(x(floor(N_x/2)+1:end),abs(psi_store(floor(N_x/2)+1:end,end)).^2);
        Pop_distribution_t(:,BFGS_iter)     = [N1_t;N2_t];
        a_opt_vals(BFGS_iter)               = a_opt;

        % if after first step, check if step was not effective
        if BFGS_iter > 1
            if abs((cost_vals(BFGS_iter) - cost_vals(BFGS_iter - 1))/cost_vals(BFGS_iter - 1)) < 5e-4
                NE_steps = NE_steps + 1;
                fprintf('reduction factor = %e\n',abs((cost_vals(BFGS_iter) - cost_vals(BFGS_iter - 1))/cost_vals(BFGS_iter - 1)))
            else
                NE_steps = 0;
            end
        end
        fprintf('Step %i done in %.2f s, with lambda cost: %.2f percent, #NE = %i\n',BFGS_iter,iter_time,lambda_cost_val(BFGS_iter)/val_opt*100,NE_steps)


        % check if not effective steps reached max_NE_steps
        if NE_steps >= max_NE_steps
            breaker = 1;
        end

        % break if BFGS got ineffective
        if breaker
            cost_vals(BFGS_iter:end)        = val_opt;
            lambda_cost_val(BFGS_iter:end)  = lambda_cost(t,A_t_is);
            fprintf('too many ineffective steps\n')
            break;
        end



    end
    % plot result of current ramp
    figure(12)
    plot(t,A_t_is)
    grid on
    hold on
    xlabel('t in ms')
    ylabel('A')

    fprintf('cost function reduction factor: %f\n',initial_cost/cost_vals(end))

    % store combined ramp if ramp was used to create initial state
    if switch_scenario == 4
        A_t_is_cases{indx_ramp}                 = [initial_state_struct.A_t_previous(1:end-1),A_t_is];
        t_case{indx_ramp}                       = [initial_state_struct.t(1:end-1),initial_state_struct.t(end)+t_case_ramp{indx_ramp}];
    else
        A_t_is_cases{indx_ramp}                 = A_t_is;
        t_case{indx_ramp}                       = t_case_ramp{indx_ramp};
    end
    
    % lambda_cost_cases(:,indx_ramp)          = lambda_cost_val;
    % cost_vals_cases(:,indx_ramp)            = [initial_cost,cost_vals];
    % iter_times_cases(indx_ramp,:)           = iter_times;
    % iter_n_iters_opt_cases(indx_ramp,:)     = iter_n_iters_opt;
    % a_opt_vals_cases(:,indx_ramp)           = a_opt_vals;

end
