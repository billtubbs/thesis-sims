%% Run observer simulations
%
% Simulate different observers on measurement data from grinding
% simulation model
%
% Author: Bill Tubbs
%
% Input files:
%  - rod_obs_P1D_c4.m - process model and observers
%
% Input data:
%  - Simulation outputs from Simulink model simulations.
%    See data folder 
%
% Results files:
%  1. rod_obs_sim_1_2.csv - time series data
%  2. rod_obs_sim_1_2_MMKF.csv - additional time series data
%  3. rod_obs_sim_1_summary.csv - simulation parameters
%  4. rod_obs_sim_1_resps.csv - observer responses to disturbances
%
% Notes:
%  1. When running this file some existing results are over-
%     written, whereas others are appended to.  To start with a
%     'blank slate', delete all the files in the 'results' folder
%     before running this script.
%
%     The files rod_obs_sim_1_summary.csv and rod_obs_sim_1_resps.csv
%     contain accumulated results of previous simulations.
%

clear all

% Specify path to observer functions and others
addpath('~/process-observers')
addpath('../data-utils')
addpath('../plot-utils')

% Sub-directories used
data_dir = 'data';
results_dir = 'results/sens';
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% Specify application case
%p_case = 1;  % Vary process model parameters
p_case = 2;  % Vary RODD model parameters
 
% Specify which data set to run simulations with
% I used:
% - 1 for process model estimation (Fig. 4 in paper)
% - 2 for process model validation (model selection)
% - 3 for initial observer test (Fig. 5 in paper)
% - 5 for observer parameter optimization
% - 6 to 15 for observer Monte Carlo simulations.
% - 6 for sensitivity analyses.
i_in_seq = 6;

switch p_case

    case 1

        % Range of relative uncertainty of parameters
        var_ratios = [1/2 1/1.5 1/1.15 1 1.15 1.5 2];
        n_ratios = numel(var_ratios);

        % Set up grid of all param ratio combinations
        [Kp_ratios, Tp1_ratios] = ndgrid(1:n_ratios, 1:n_ratios);
        combs = [reshape(Kp_ratios, [], 1) reshape(Tp1_ratios, [], 1)];

    case 2
        
        % Range of relative uncertainty of parameters
        var_ratios = [1/2 1/1.5 1/1.15 1 1.15 1.5 2];
        n_ratios = numel(var_ratios);

        % Set up grid of all param ratio combinations
        [epsilon_ratios, sigma_ratios] = ndgrid(1:n_ratios, 1:n_ratios);
        combs = [reshape(epsilon_ratios, [], 1) reshape(sigma_ratios, [], 1)];

end

% Run observer simulations
n_combs = size(combs, 1);
for i_iter = 1:n_combs

    % Combination of variations for this iteration
    comb = combs(i_iter, :);

    fprintf("\nStarting observer simulation #%d of %d...\n", i_iter, n_combs)

    % Specify observer model parameters
    switch p_case

        case 1

            % Kp : default -32.4
            % Tp1 : default 0.106
            % Tp2 : default 0.106
            % thetap : default 0.2
            Kp = -32.4 * var_ratios(comb(1));  % adjusted
            Tp1 = 0.106 * var_ratios(comb(2));  % adjusted
            Tp2 = Tp1;  % same as Tp1 (adjusted)
            thetap = 0.2;  % fixed
            epsilon = 0.01;  % default
            sigma_wp_1 = 0.2717 / 100;  % default
            b = 100;  % default

        case 2

            Kp = -32.4;  % default
            Tp1 = 0.106;  % default
            Tp2 = 0.106;  % default
            thetap = 0.2;  % default
            % epsilon : default 0.01
            % step_mag : default 0.2717
            % b : default 100
            % sigma_wp_1 : default 0.2717 / 100;
            % sigma_wp_2 : default 0.2717
            epsilon = 0.01 * var_ratios(comb(1));  % adjusted
            sigma_wp_1 = 0.2717 / 100;
            b = 100 * var_ratios(comb(2));  % adjusted

    end

    % Load observers
    % These are adjusted according to parameters specified above

    rod_obs_P2DcTd4_adj  % model identified from data
    %rod_obs_P2Dcd1_T_adj  % ident. from true outputs

    % Choose which observers to simulate
    % - KF1 : Kalman filter tuned to minimize steady-state errors
    % - KF3 : Kalman filter tuned to minimize overall errors
    % - SKF : Scheduled Kalman filter
    % - MKF_SF95 : Multiple model observer - sequence fusion
    % - MKF_SF1 : Multiple model observer - sequence fusion
    % - MKF_SP1 : Multiple model observer - sequence pruning
    % - SKF : Scheduled Kalman filter

    observers = {KF1, KF3, MKF_SF95, MKF_SF1, MKF_SP1, SKF};

    % Load system simulation results
    if i_in_seq < 6
        nT = 300;
    else
        nT = 2460;
    end
    filename = sprintf('sim_OL_rc_est_mix_factor_%d_%d_ident.csv', ...
                       nT, i_in_seq);
    sim_data = readtable(fullfile(data_dir, filename));

    t = sim_data.t;
    t_stop = t(end);
    nT = ceil(t_stop / Ts);
    assert(size(t, 1) == nT+1)

    U = zeros(nT+1, 0);
    Pd = sim_data{:, 'BASE_ORE_MIX'};
    Y = sim_data{:, 'SAG_OF_P80'};
    Y_m = sim_data{:, 'SAG_OF_P80_M'};  % with measurement noise

    % Calculate random shock signal that would replicate the
    % simulated disturbance
    n_dist = size(Pd, 2);
    Wp = [diff(Pd); zeros(1, n_dist)];  % shifted for delay
    assert(isequal(size(Wp), [nT+1 n_dist]))
    % Find when shocks occurred TODO: should generate these at the
    % time the simulations are run.
    [rows,cols,v] = find(Wp);
    alpha = zeros(nT+1, n_dist);
    for i = 1:numel(rows)
        alpha(rows(i), cols(i)) = 1;
    end
    if n_dist == 1
        gamma = alpha;
    end
    assert(n_dist == 1)

    % Calculate plant output predictions with the model
    Y_model = lsim(Gpss, Wp, t);

    % To test observers, try the following
    % (i) remove measurement noise
    %Y_m = Y;
    % (ii) set Y_m to model predictions
    %Y = Y_model;
    %Y_m = Y;

    % Run simulation
    input_data = table(U, alpha, gamma, Wp, Pd, Y, Y_m);
    sim_out = run_obs_simulation(Ts, input_data, observers);
    observers = sim_out.observers;  % Updated observers


    %% Display and save simulation results
    
    % Remove semi-colon to display results table
    sim_out.data

    filename = sprintf('rod_obs_sim_%d_%d.csv', p_case, i_in_seq);
    writetable(drop_empty_cols(sim_out.data), fullfile(results_dir, filename));
    fprintf("Observer simulation results saved to file: %s\n", filename)

    % Count number of observers and MMKF observers
    n_obs = numel(observers);
    n_obs_mkf = 0;
    observers_mkf = double.empty(1, 0);
    for i=1:n_obs
        if startsWith(observers{i}.label, "MMKF")
            n_obs_mkf = n_obs_mkf + 1;
            observers_mkf(n_obs_mkf) = i;
        end
    end

    % Simulation results into variables
    t = sim_out.data{:,'t'};
    U = sim_out.data{:,'U'};
    alpha = sim_out.data{:, 'alpha'};
    X_est = sim_out.data{:, vector_element_labels('X_est', '', n_obs)};
    Y = sim_out.data{:, 'Y'};
    Y_m = sim_out.data{:, 'Y_m'};
    Y_est = sim_out.data{:, vector_element_labels('Y_est', '', n_obs)};
    E_obs = sim_out.data{:, 'E_obs'};

    % Save results from multiple model filters (if used)
    for f = 1:n_obs_mkf

        label = observers{observers_mkf(f)}.label;
        MKF_sim_results = [sim_out.data(:, {'k', 't'}) ...
            array2table_with_name(sim_out.MKF_i{f}, 'i', '_') ...
            array2table_with_name(sim_out.MKF_p_seq_g_Yk{f}, 'p_seq_g_Yk', '_') ...
            array2table_with_name(sim_out.MKF_X_est{f}, 'X_est', '_') ...
        ];

        filename = sprintf('rod_obs_sim_%d_%d_%s.csv', p_case, i_in_seq, label);
        writetable(drop_empty_cols(MKF_sim_results), fullfile(results_dir, filename));

    end
    
    fprintf("MKF simulation results saved to file: %s\n", filename)


    %% Prepare labels for tables and plots

    rod_obs_make_labels


    %% Compute observer performance metrics

    % The following script uses the values stored in
    [metrics, metrics_params, errors, metrics_labels] = calculate_obs_metrics(Y, Y_est, ...
        obs_labels, Pd, Ts);

    % Make metrics labels for all observers, e.g. for observer 'KF1':
    %  - 'MSE_y_est_KF1' : overall MSE
    %  - 'MSE_tr_y_est_KF1' : MSE in transition periods
    %  - 'MSE_ss_y_est_KF1' :  MSE in steady-state periods
    %  - 'Var_ss_y_est_KF1' : Variance in steady-state periods
    n_metrics = numel(metrics_labels);
    obs_metrics_labels = cell(n_metrics, n_obs * ny);
    for i = 1:n_metrics
        metric_label = metrics_labels{i};
        labels = matrix_element_labels(metric_label, y_est_labels, obs_labels, '');
        obs_metrics_labels(i, :) = labels(:)';
    end

    %% Display RMSE results

    % Transpose the table (complicated in MATLAB):
    rmse_table_tr = rows2vars(metrics);
    rmse_table_tr = removevars(rmse_table_tr, 'OriginalVariableNames');
    rmse_table_tr.Properties.RowNames = {'RMSE', ...
        'RMSE in transitions', 'RMSE in steady-state', ...
        'Variance in steady-state', 'RMSD in steady-state'};
    disp(rmse_table_tr)

    % Compute errors in MKF observer estimates (if used)
    MKF_Y_errors = cell(size(sim_out.MKF_Y_est));
    MKF_Y_RMSE = cell(1, n_obs_mkf);
    for f = 1:n_obs_mkf
        obs = observers{observers_mkf(f)};
        MKF_Y_RMSE{f} = size(sim_out.MKF_Y_est{f}, 1);
        % Compute errors in multiple filter state estimates
        % Find out how many hypotheses were saved
        nh = size(sim_out.MKF_p_seq_g_Yk{f}, 2);
        MKF_Y_errors{f} = repmat(Y, 1, nh) - sim_out.MKF_Y_est{f};
        MKF_Y_RMSE{f} = mean(MKF_Y_errors{f}.^2, 1);
    end


    %% Combine all parameters and results and add to summary results file

    % System model parameters
    sys_params = objects2tablerow(containers.Map({'sys'}, {model}));

    % Observer parameters
    rv_params = objects2tablerow( ...
        containers.Map({'epsilon', 'sigma_wp', 'sigma_M'}, ...
            {epsilon, sigma_wp, sigma_M}) ...
    );
    obs_params = cell(1, n_obs);
    for f = 1:n_obs
        obs = observers{f};
        params = get_obs_params(obs);
        objects = containers.Map(cellstr(obs.label), {params});
        obs_params{f} = objects2tablerow(objects);
    end
    obs_params = horzcat(obs_params{:});

    % Simulation settings
    sim_params = table(p_case, i_in_seq, t_stop, Ts, nT, nu, ny, n_obs);

    % Observer metrics
    obs_metrics = [ ...
        objects2tablerow(containers.Map({'metrics'}, {metrics_params})) ...
        array2table(reshape(metrics.Variables', [], ...
            n_obs*n_metrics), 'VariableNames', obs_metrics_labels);
    ];

    % Add additional variables from process model
    if p_case == 1
        extra_vars = table(Kp, Tp1);
    else
        extra_vars = table();
    end

    % Summary table
    summary_results = [ ...
        array2tablerow(datetime(), 'Time') ...
        extra_vars ...
        sim_params ...
        sys_params ...
        array2tablerow(obs_labels, 'obs') ...
        rv_params ...
        obs_params ...
        obs_metrics ...
    ];

    % Save to csv file
    filename = sprintf('rod_obs_sim_sens_%d_summary.csv', p_case);
    if isfile(fullfile(results_dir, filename))
        % Load existing results and combine
        resp_data_existing = readtable(fullfile(results_dir, filename));
        fprintf("Existing results loaded from file: %s\n", filename)
        summary_results = vertcat(resp_data_existing, summary_results);
    end

    % Save all results to file
    writetable(drop_empty_cols(summary_results), fullfile(results_dir, filename));
    fprintf("Summary results saved to file: %s\n", filename)


     %% Compute disturbance response trajectories
%     % remove to save time
% 
%     % Choose length of responses in sample periods
%     nT_resp = 51;  % 0 to 1.5 hrs
% 
%     % Find step changes in disturbance
%     [idxs, diffs] = find_step_periods(Pd, metrics_params.n_settle, nT_resp);
%     n_resp = numel(idxs);
%     fprintf("Step responses identified: %d\n", n_resp)
%     Y_resps = nan(nT_resp, n_resp);
%     Y_est_resps = repmat({nan(nT_resp, n_resp*1)}, 1, n_obs);
%     for i = 1:n_resp
%         idx = idxs{i};
%         nT_max = min(diff(idx)+1, nT_resp);
%         Y_resps(:, i) = Y(idx(1):idx(1)+nT_max-1);
%         for f = 1:n_obs
%             Y_est_resps{f}(:, i) = Y_est(idx(1):idx(1)+nT_max-1, ...
%                 (f-1)*ny+1:f*ny);
%         end
%     end
%     % Convert Y_est_resp_obs to tables
%     for f = 1:n_obs
%         label = strcat('yk_est_', obs_labels{f});
%         Y_est_resps{f} = array2table_with_name(Y_est_resps{f}', label, '_');
%     end
%     resp_data = [ ...
%         array2table_with_name(repmat(i_in_seq, n_resp, 1), 'i_in_seq', '_') ...
%         table(diffs) ...
%         array2table_with_name(Y_resps', 'yk', '_') ...
%         horzcat(Y_est_resps{:}) ...
%     ];
% 
%     % Save to csv file
%     filename = sprintf('rod_obs_sim_%d_resps.csv', p_case);
%     if isfile(fullfile(results_dir, filename))
%         % Load existing results and combine
%         resp_data_existing = readtable(fullfile(results_dir, filename));
%         fprintf("Existing step responses loaded from file: %s\n", filename)
%         resp_data = vertcat(resp_data_existing, resp_data);
%     end
% 
%     % Save all results to file
%     writetable(drop_empty_cols(resp_data), fullfile(results_dir, filename));
%     fprintf("Step responses saved to file: %s\n", filename)

end


% For results plots, run this file next:
%rod_obs_sim_plots

fprintf("Run rod_obs_sim_plots.m to produce plots.\n")

% To calculate evaluation metrics, run this file:
%rod_obs_calc_metrics

fprintf("Run rod_obs_calc_metrics.m to calculate evaluation metrics.\n")
