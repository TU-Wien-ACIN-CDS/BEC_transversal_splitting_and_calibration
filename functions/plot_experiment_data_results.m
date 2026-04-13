function plot_experiment_data_results(GPE_results_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots the experiment results
% 
% Args:
%     GPE_results_struct (struct):          structure including data from experiment simulation
%
% Returns:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_ramp     = GPE_results_struct.N_ramps;

%% Ramps
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    plot(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.Szenarios{indx_ramp}.A_t)
    grid on
    xlabel('t in ms')
    ylabel('A')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'--')
end
title(t1,'Ramps')
set(gcf,'Name','Ramps')

%% density
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    imagesc(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.x,abs(GPE_results_struct.Szenarios{indx_ramp}.Psi_x_t).^2)
    xlabel('t in ms')
    ylabel('x in um')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'w--')
end
title(t1,'Density')
set(gcf,'Name','Density')

%% phase
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    imagesc(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.x,angle(GPE_results_struct.Szenarios{indx_ramp}.Psi_x_t).^2)
    xlabel('t in ms')
    ylabel('x in um')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'w--')
end
title(t1,'Phase')
set(gcf,'Name','Phase')

%% density of Momentum
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    imagesc(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.Szenarios{indx_ramp}.k,abs(GPE_results_struct.Szenarios{indx_ramp}.Psi_hat_k_t).^2)
    xlabel('t in ms')
    ylabel('x in um')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'w--')
end
title(t1,'Density of Momentum')
set(gcf,'Name','Density of Momentum')

%% phase of Momentum
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    imagesc(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.Szenarios{indx_ramp}.k,angle(GPE_results_struct.Szenarios{indx_ramp}.Psi_hat_k_t).^2)
    xlabel('t in ms')
    ylabel('x in um')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'w--')
end
title(t1,'Phase of Momentum')
set(gcf,'Name','Phase of Momentum')


%% Interwell distance
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    plot(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.Szenarios{indx_ramp}.d_t)
    grid on
    xlabel('t in ms')
    ylabel('d in um')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'--')
end
title(t1,'Interwell distance')
set(gcf,'Name','Interwell distance')


%% Sigma in situ
figure
t1 = tiledlayout('flow');
for indx_ramp = 1:N_ramp
    nexttile
    plot(GPE_results_struct.Szenarios{indx_ramp}.t,GPE_results_struct.Szenarios{indx_ramp}.sigma_t)
    grid on
    xlabel('t in ms')
    ylabel('sigma in um')
    xline(GPE_results_struct.Szenarios{indx_ramp}.T_ramp,'--')
end
title(t1,'Sigma in situ')
set(gcf,'Name','Sigma in situ')

%%

d_amp_list      = zeros(1,N_ramp);
d_freq_list        = zeros(1,N_ramp);
d_mean_list     = zeros(1,N_ramp);
d_var_list      = zeros(1,N_ramp);

sigma_amp_list      = zeros(1,N_ramp);
sigma_freq_list        = zeros(1,N_ramp);
sigma_mean_list     = zeros(1,N_ramp);
sigma_var_list      = zeros(1,N_ramp);

delta_E_list        = zeros(1,N_ramp);

for indx_ramp = 1:N_ramp
    d_amp_list(indx_ramp)               = GPE_results_struct.Szenarios{indx_ramp}.d_amp;
    d_freq_list(indx_ramp)              = GPE_results_struct.Szenarios{indx_ramp}.d_freq;
    d_mean_list(indx_ramp)              = GPE_results_struct.Szenarios{indx_ramp}.d_mean;
    d_var_list(indx_ramp)               = GPE_results_struct.Szenarios{indx_ramp}.d_var;

    sigma_amp_list(indx_ramp)           = GPE_results_struct.Szenarios{indx_ramp}.sigma_amp;
    sigma_freq_list(indx_ramp)          = GPE_results_struct.Szenarios{indx_ramp}.sigma_freq;
    sigma_mean_list(indx_ramp)          = GPE_results_struct.Szenarios{indx_ramp}.sigma_mean;
    sigma_var_list(indx_ramp)           = GPE_results_struct.Szenarios{indx_ramp}.sigma_var;

    delta_E_list(indx_ramp)             = GPE_results_struct.Szenarios{indx_ramp}.delta_E;
end

% interwell distance results
figure
subplot(5,1,1)
plot(1:N_ramp,d_amp_list,'*--')
xlabel('Ramp Number')
ylabel('amp of d in um')
title('Interwell distance values: amp')
xlim([0,N_ramp+1])
grid on
subplot(5,1,2)
plot(1:N_ramp,d_freq_list,'*--')
xlabel('Ramp Number')
ylabel('frequency in kHz')
title('frequency')
xlim([0,N_ramp+1])
grid on
subplot(5,1,3)
plot(1:N_ramp,d_mean_list,'*--')
xlabel('Ramp Number')
ylabel('mean d in um')
title('mean')
xlim([0,N_ramp+1])
grid on
subplot(5,1,4)
plot(1:N_ramp,d_var_list,'*--')
xlabel('Ramp Number')
ylabel('var d in um^2')
title('Var')
xlim([0,N_ramp+1])
grid on
subplot(5,1,5)
plot(1:N_ramp,delta_E_list,'*--')
xlabel('Ramp Number')
ylabel('energy over groundstate')
title('delta Energy')
xlim([0,N_ramp+1])
grid on
set(gcf,'Name','Interwell distance values')

% sigma in situ results
figure
subplot(5,1,1)
plot(1:N_ramp,sigma_amp_list,'*--')
xlabel('Ramp Number')
ylabel('amp of sigma in um')
title('Sigma values: amp')
xlim([0,N_ramp+1])
grid on
subplot(5,1,2)
plot(1:N_ramp,sigma_freq_list,'*--')
xlabel('Ramp Number')
ylabel('frequency in kHz')
title('frequency')
xlim([0,N_ramp+1])
grid on
subplot(5,1,3)
plot(1:N_ramp,sigma_mean_list,'*--')
xlabel('Ramp Number')
ylabel('mean sigma in um')
title('mean')
xlim([0,N_ramp+1])
grid on
subplot(5,1,4)
plot(1:N_ramp,sigma_var_list,'*--')
xlabel('Ramp Number')
ylabel('var sigma in um^2')
title('Var')
xlim([0,N_ramp+1])
grid on
subplot(5,1,5)
plot(1:N_ramp,delta_E_list,'*--')
xlabel('Ramp Number')
ylabel('energy over groundstate')
xlim([0,N_ramp+1])
grid on
set(gcf,'Name','Sigma in situ values')
title('delta Energy')


end