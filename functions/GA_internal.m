%% This is an internal function, be very carefull when changing anything!!!!

%% Collecting sensitivity matrices and measurement value vectors

S_all_f         = zeros(N_all,N_par);        % Sensitivity matrix single linear ramp frequencies pre norm
S_all_amp       = zeros(N_all,N_par);        % Sensitivity matrix single linear ramp amplitudes pre norm
S_all_d         = zeros(N_all,N_par);        % Sensitivity matrix single linear ramp mean differences pre norm
S_all_ph        = zeros(N_all,N_par);        % Sensitivity matrix single linear ramp phase pre norm
S_all_off       = sensitivity_offset_all;    % Sensitivity matrix offset sloshing frequencies pre norm
S_f_ML          = zeros(N_ML,N_par);         % Sensitivity matrix multiple linear ramp frequencies pre norm
S_amp_ML        = zeros(N_ML,N_par);         % Sensitivity matrix multiple linear ramp amplitudes pre norm
S_d_ML          = zeros(N_ML,N_par);         % Sensitivity matrix multiple linear ramp mean difference pre norm
S_ph_ML         = zeros(N_ML,N_par);         % Sensitivity matrix multiple linear ramp phases pre norm
f_all           = zeros(N_all,1);            % frequencies single linear ramps
amp_all         = zeros(N_all,1);            % amplitudes single linear ramps
A_end_all_lin   = zeros(N_all,1);            % end parameters single linear ramps
T_all_lin       = zeros(N_all,1);            % ramp times single linear ramps
f_ML_all        = zeros(N_ML,1);             % frequencies multiple linear ramps
amp_ML_all      = zeros(N_ML,1);             % amplitudes multiple linear ramps
ph_ML_all       = zeros(N_ML,1);             % phases multiple linear ramps
A_end_ML_all    = zeros(N_ML,1);             % end parameters multiple linear ramps
T_end_ML_all    = zeros(N_ML,1);             % combined ramp times multiple linear ramps

norm_ver = 1; % selector for which norm is used

% collecting single linear ramp values
for indx_T  = 1:N_T
    for indx_A = 1:N_A
        indx_all                = N_A*(indx_T-1) + indx_A;
        S_all_f(indx_all,:)     = sensitivity_f_all(indx_A,indx_T,:);
        S_all_amp(indx_all,:)   = sensitivity_amp_all(indx_A,indx_T,:);
        S_all_d(indx_all,:)     = sensitivity_d_all(indx_A,indx_T,:);
        S_all_ph(indx_all,:)    = sensitivity_ph_all(indx_A,indx_T,:);
        f_all(indx_all)         = f_grid_all(indx_A,indx_T);
        amp_all(indx_all)       = amp_grid_all(indx_A,indx_T);
        A_end_all_lin(indx_all)     = A_end_grid_all(indx_A,indx_T);
        T_all_lin(indx_all)         = T_grid_all(indx_A,indx_T);
    end
end

% collecting multiple linear ramp values
for indx_A_ML = 1:N_A_ML
    for indx_T_ML = 1:N_T_ML
        indx_all_ML                 = (indx_A_ML-1)*N_T_ML + indx_T_ML;
        S_f_ML(indx_all_ML,:)       = sensitivity_f_ML_all(indx_A_ML,indx_T_ML,:);
        S_amp_ML(indx_all_ML,:)     = sensitivity_amp_ML_all(indx_A_ML,indx_T_ML,:);
        S_ph_ML(indx_all_ML,:)      = sensitivity_ph_ML_all(indx_A_ML,indx_T_ML,:);
        f_ML_all(indx_all_ML,:)     = f_ML(indx_A_ML,indx_T_ML);
        amp_ML_all(indx_all_ML,:)   = amp_ML(indx_A_ML,indx_T_ML);
        ph_ML_all(indx_all_ML,:)    = ph_ML(indx_A_ML,indx_T_ML);
        A_end_ML_all(indx_all_ML,:) = A_end_grid_ML(indx_A_ML,indx_T_ML);
        T_end_ML_all(indx_all_ML,:) = T_grid_ML(indx_A_ML,indx_T_ML);
    end
end

S_cells                 = {S_all_f,sensitivity_offset_all,S_f_ML};                        % cell with all loaded sensitivity
A_all_cells             = {A_end_all_lin,A_offset_list,A_end_ML_all};          % cell with all loaded A values
T_all_cells             = {T_all_lin,zeros(size(A_offset_list)),T_end_ML_all};     % cell with all loaded T values
type_selection          = [1,2,3];                                                                                  % type of measurement for sensitivity
type_names              = {'Frequency of single linear ramp','Frequency of offset sloshing','Frequency of multiple linear ramps'};
type_letters            = {'f','f','f'};                                                                % letter symbalizing the measurements
measurement_types       = {'SL','O','ML'};                                                        % underlying measurement method
N_type_measurements     = numel(S_cells);                                                                                       % number of sensitivities
S_all                   = [];                                                                                                   % collected Sensititivity matrix
running_indx            = 0;                                                                                                    % running index for correct selector definition
selector_cells          = cell(N_type_measurements,1);                                                                          % selectors for different sensitivities
vecnorm_cells           = cell(N_type_measurements,1);                                                                          % vectornorm for different sensitivities

par_names       = {'k1','k2','c','As','g1D'};

T_all                   = [];
A_all                   = [];

for indx_cells  = 1:N_type_measurements
        [N_current,~]               = size(S_cells{indx_cells});
        switch norm_ver
            case 0  % no normalization
                S_current                   = S_cells{indx_cells};
                vecnorm_cells{indx_cells}   = 1;
                S_current_norm              = S_current ./ vecnorm_cells{indx_cells};
            case 1  % divided by vectornorm of collums
                S_current                   = S_cells{indx_cells};
                vecnorm_cells{indx_cells}   = vecnorm(S_current,2,1);
                S_current_norm              = S_current ./ vecnorm_cells{indx_cells};
            case 2 % divided by absolute maximum of collums
                S_current                   = S_cells{indx_cells};
                vecnorm_cells{indx_cells}   = max(abs(S_current),[],'all');
                S_current_norm              = S_current ./ vecnorm_cells{indx_cells};
            case 3 % normalized to [0 1]
                S_current                   = S_cells{indx_cells};
                S_current_norm              = S_current - min(S_current,[],'all');
                vecnorm_cells{indx_cells}   = max(abs(S_current_norm),[],'all');
                S_current_norm              = S_current_norm ./ vecnorm_cells{indx_cells};
            case 4  % no normalization
                S_current                   = S_cells{indx_cells};
                vecnorm_cells{indx_cells}   = norm_vals_measurements(indx_cells);
                S_current_norm              = S_current ./ vecnorm_cells{indx_cells};
        end
        S_all                       = [S_all;S_current_norm];
        selector_cells{indx_cells}  = running_indx + (1:N_current);
        running_indx                = running_indx + N_current;
        T_all       = [T_all;T_all_cells{indx_cells}(:)];
        A_all       = [A_all;A_all_cells{indx_cells}(:)];
end

N_used      = numel(A_all);

%% plotting selected sensitivities
indx_lims = [-inf inf];



nexttile
plot(A_all)
grid on
xlim(indx_lims)
ylabel('A end')
nexttile
plot(T_all)
grid on
xlim(indx_lims)
ylabel('T')

nexttile
contourf(T_list,A_end_list,f_grid_all)
xlabel('T in ms')
ylabel('A')
title('f')
nexttile
contourf(T_list,A_end_list,amp_grid_all)
xlabel('T in ms')
ylabel('A')
title('amp')

%% implement GA here

overall_best        = 1;
overall_best_cost   = inf;

N_start             = 30;
N_population        = 240;

p_T                 = 0.001;
p_plus              = 0.001;
p_minus             = 0.01;

iter_GA_max         = 1e4;

population          = cell(1,N_population);
for indx_pop = 1:N_population
    pop                         = unique(randi([1,N_used],N_start,1));
    N_current                   = numel(pop);
    population{indx_pop}        = pop(randperm(N_current,N_current));
end

min_val_GA          = zeros(1,iter_GA_max);
best_val_GA         = min_val_GA;

times_store         = zeros(1,iter_GA_max);
mean_pop_store      = zeros(1,iter_GA_max);

for iter_GA = 1:iter_GA_max
    tic
    % calculating fitness
    fitness    = zeros(1,N_population);
    size_pops   = zeros(N_population,1);
    for indx_pop = 1:N_population
        S_pop               = S_all(population{indx_pop},:);
        fitness(indx_pop)   = calculate_rho_sens_from_sensitivity_matrix(S_pop);
        size_pops(indx_pop) = numel(population{indx_pop});
    end
    mean_pop_store(iter_GA) = mean(size_pops);

    % checking for new lowest cost/fitness
    [min_val,min_indx]      = min(fitness);
    if min_val < overall_best_cost
        overall_best        = population{min_indx};
        overall_best_cost   = min_val;

        fprintf('iteration %i best is: %.3f\n',iter_GA,overall_best_cost)
        % if min_val < 3
        %     break;
        % end
    end
    min_val_GA(iter_GA:end)     = min_val;
    best_val_GA(iter_GA:end)    = overall_best_cost;

    % THE CULLING (population reduction)
    [~,indx_sorted]     = sort(fitness);
    indx_selected       = zeros(size(fitness));
    switch CULLING_Variant
        case 1
            indx_selected(indx_sorted(1:N_population/3)) = 1;
        case 2
            indx_selected(indx_sorted(randperm(N_population,N_population/3))) = 1;
        case 3
            p_min                   = 0.3;
            p_max                   = 0.95;
            p_X                     = linspace(p_min,p_max,N_population);
            rand_X                  = rand(size(p_X));
            diff_p                  = p_X - rand_X;
            [~,indx_diff_sorted]    = sort(diff_p);
            indx_selected(indx_sorted(indx_diff_sorted(1:N_population/3))) = 1;
    end
    indx_notselected    = ~indx_selected;

    places_selected     = find(indx_selected);
    places_not_selected = find(indx_notselected);

    places_selected     = places_selected(randperm(N_population/3,N_population/3));
    places_not_selected = places_not_selected(randperm(2*N_population/3,2*N_population/3));
    
    % crossover)
    for indx_culling = 1:numel(places_selected)
        parent1 = population{places_selected(indx_culling)};
        parent2 = population{end-places_selected(indx_culling)+1};

        N1      = numel(parent1);
        N2      = numel(parent2);

        collected_parents   = union(parent1,parent2);
        N_parents           = numel(collected_parents);

        Nc1     = round((N1+N2)/2+2*(rand(1)-0.5)*abs(N1-N2)/2);
        if Nc1 > N_parents
            child1  = [collected_parents(randperm(N_parents,N_parents));randperm(N_used,Nc1-N_parents)'];
        else
            child1  = collected_parents(randperm(N_parents,Nc1));
        end
        Nc2     = round((N1+N2)/2+2*(rand(1)-0.5)*abs(N1-N2)/2);
        if Nc2 > N_parents
            child2  = [collected_parents(randperm(N_parents,N_parents));randperm(N_used,Nc2-N_parents)'];
        else
            child2  = collected_parents(randperm(N_parents,Nc2));
        end

        population{places_not_selected(indx_culling)*2-1}   = child1;
        population{places_not_selected(indx_culling)*2}     = child2;

    end

    % Mutation
    switch mutation_var
        case 1 % mutation of all members
            for indx_pop = 1:N_population
                current_pop = population{indx_pop};
                if isempty(current_pop)
                    current_pop = randperm(N_used,1);
                end
                for indx_element = 1:numel(current_pop)
                    if rand(1) < p_T
                        current_pop(indx_element) = randperm(N_used,1); 
                    end
                end
                if rand(1) < p_minus
                    current_pop(randperm(numel(current_pop),1)) = [];
                end
                if rand(1) < p_plus
                    current_pop = [current_pop;randperm(N_used,1)];
                end
                population{indx_pop}    = unique(current_pop);

            end
        case 2 % mutation only of new members
            for indx_new = 1:numel(places_not_selected)
                current_pop = population{places_not_selected(indx_new)};
                if isempty(current_pop)
                    current_pop = randperm(N_used,1);
                end
                for indx_element = 1:numel(current_pop)
                    if rand(1) < p_T
                        current_pop(indx_element) = randperm(N_used,1);
                    end
                end
                if rand(1) < p_minus
                    if numel(current_pop) > 1
                        current_pop(randperm(numel(current_pop),1)) = [];
                    else
                        current_pop = randi(N_used,1);
                    end
                end
                if rand(1) < p_plus
                    current_pop = [current_pop;randperm(N_used,1)];
                end
                population{places_not_selected(indx_new)}   = unique(current_pop);
            end
    end

    if mod(iter_GA,iter_GA_max/10) == 0
        fprintf('%.0f percent done\n',iter_GA/iter_GA_max*100)
    end
    times_store(iter_GA) = toc;
end

% figure
% subplot(2,1,1)
% semilogy(min_val_GA)
% grid on
% hold on
% semilogy(best_val_GA,'k--')
% subplot(2,1,2)
% plot(mean_pop_store)
% grid on


selected_experiments    = overall_best;

disp(mean(times_store))

%%
N_selected              = numel(selected_experiments);

indx1_selected          = zeros(N_type_measurements,N_selected);
indx2_selected          = zeros(N_type_measurements,N_selected);
allsizes_minus          = zeros(N_type_measurements,1);
indx1_cell              = cell(N_type_measurements,1);
indx2_cell              = cell(N_type_measurements,1);

for indx_cell = 1:N_type_measurements
    allsizes_minus(indx_cell)    = numel(selector_cells{indx_cell});
end
current_size            = cumsum([0,;allsizes_minus]);

for indx_exp = 1:N_selected
    indx_all_min = selected_experiments(indx_exp);
    for indx_cell = 1:N_type_measurements
        if ismember(indx_all_min,selector_cells{indx_cell})
            switch type_selection(indx_cell)
                case 1
                    [indx1,indx2]                       = ind2sub([N_A,N_T],indx_all_min - current_size(indx_cell));
                    indx1_selected(indx_cell,indx_exp)  = indx1;
                    indx2_selected(indx_cell,indx_exp)  = indx2;
                    indx1_cell{indx_cell}               = [indx1_cell{indx_cell},indx1];
                    indx2_cell{indx_cell}               = [indx2_cell{indx_cell},indx2];
                case 2
                    [indx1,indx2]                       = ind2sub([N_A_offset,1],indx_all_min - current_size(indx_cell));
                    indx1_selected(indx_cell,indx_exp)  = indx1;
                    indx2_selected(indx_cell,indx_exp)  = indx2;
                    indx1_cell{indx_cell}               = [indx1_cell{indx_cell},indx1];
                    indx2_cell{indx_cell}               = [indx2_cell{indx_cell},indx2];
                case 3
                    [indx1,indx2]                       = ind2sub([N_A_ML,N_T_ML],indx_all_min - current_size(indx_cell));
                    indx1_selected(indx_cell,indx_exp)  = indx1;
                    indx2_selected(indx_cell,indx_exp)  = indx2;
                    indx1_cell{indx_cell}               = [indx1_cell{indx_cell},indx1];
                    indx2_cell{indx_cell}               = [indx2_cell{indx_cell},indx2];
            end
        end
    end


end

%%
S_selected      = S_all(selected_experiments,:);
S_w_selected    = (S_selected ./ vecnorm(S_selected));

S_w_selected_square = S_w_selected'*S_w_selected;

I_CL_selected   = 1./sqrt(min(abs(eig(S_w_selected_square))));
I_Sw_selected   = cond(S_w_selected_square);

index_text      = sprintf('I_{CL} = %.3f, I_{Sw} = %.3f, N_{exp} = %i\n',I_CL_selected,I_Sw_selected,N_selected);
fprintf(index_text)

[A_end_list_mesh,T_list_mesh]       = meshgrid(A_end_list,T_list);
[A_end_list_ML_mesh,T_end_list_ML_mesh]       = meshgrid(A_end_list_ML,T_end_list_ML);

% figure
% tl = tiledlayout("flow","TileSpacing","tight");
% for indx_cell = 1:N_type_measurements
%     if ~isempty(selector_cells{indx_cell})
%         nexttile
%         switch type_selection(indx_cell)
%             case 1
%                 scatter(A_end_list_mesh(:),T_list_mesh(:),'.','MarkerEdgeColor','#080808')
%                 hold on
%                 scatter(A_end_list(indx1_cell{indx_cell}),T_list(indx2_cell{indx_cell}),'filled','MarkerFaceColor',"#0072BD")
%                 % xlim([A_end_list(1),A_end_list(end)])
%                 % ylim([T_list(1),T_list(end)])
%             case 2
%                 scatter(A_offset_list,zeros(size(A_offset_list)),'.','MarkerEdgeColor','#080808')
%                 hold on
%                 scatter(A_offset_list(indx1_cell{indx_cell}),0*indx2_cell{indx_cell},'filled','MarkerFaceColor',"#0072BD")
%                 % xlim([A_offset_list(1),A_offset_list(end)])
%             case 3
%                 scatter(A_end_ML_all(:),T_end_list_ML(:),'.','MarkerEdgeColor','#080808')
%                 hold on
%                 scatter(A_end_ML_all(indx1_cell{indx_cell}),T_end_list_ML(indx2_cell{indx_cell}),'filled','MarkerFaceColor',"#0072BD")
%                 % xlim([A_end_ML_all(1),A_end_ML_all(end)])
%                 % % xlim([A_end_list(1),A_end_list(end)])
%         end
%         xlim([A_all(selector_cells{indx_cell}(1))-0.05,A_all(selector_cells{indx_cell}(end))+0.05])
%         ylim([T_all(selector_cells{indx_cell}(1))-0.1,T_all(selector_cells{indx_cell}(end))+0.1])
%         grid on
%         title(type_names{indx_cell})
%     end
% end
% title(tl,sprintf('I_{CL} = %.3f, I_{Sw} = %.3f',I_CL_selected,I_Sw_selected))

A_linear_ramps              = A_end_list(indx1_cell{1});
T_linear_ramps              = T_list(indx2_cell{1});
A_offset                    = A_offset_list(indx1_cell{2});
A_consecutive_linear_ramps  = A_end_list_ML(indx1_cell{3});