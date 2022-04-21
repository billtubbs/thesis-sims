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
addpath('process-observers')
addpath('~/ml-data-utils')
addpath('~/ml-plot-utils')

% Sub-directories used
data_dir = 'data';
results_dir = 'results';
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% Specify application case
p_case = 1;  % Only case 1 used here
 
% Specify which data set(s) to run simulations with
% I used:
% - 1 for process model estimation (Fig. 4 in paper)
% - 2 for process model validation (model selection)
% - 3 for initial observer test (Fig. 5 in paper)
% - 5 for observer parameter optimization
% - 6 to 15 for observer Monte Carlo simulations.
i_in_seqs = 3;
%i_in_seqs = [1, 2, 3, 4, 5];
%i_in_seqs = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

% Run observer simulations
n_in_seqs = numel(i_in_seqs);
for i_seq = 1:n_in_seqs

    i_in_seq = i_in_seqs(i_seq);
    fprintf("\nStarting observer simulations with input seq. #%d ...\n", ...
        i_in_seq)

    % Load observers
    %rod_obs_P1Dcd4
    %rod_obs_P1
    %rod_obs_P1DcD5
    %rod_obs_P2U
    %rod_obs_P2DcTd4  % observers used in draft of paper
    %rod_obs_oe125
    rod_obs_P2Dcd1_T  % ident. from true outputs

    % Choose which observers to simulate
    % - KF1 : Kalman filter tuned to minimize steady-state errors
    % - KF2 : Kalman filter tuned to minimize overall errors
    % - MMKF : Multi-model Kalman filter observer
    % - SKF : Scheduled Kalman filter

    observers = {KF1, KF2, MMKF, SKF};

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
        fprintf("MKF simulation results saved to file: %s\n", filename)

    end


    %% Prepare labels for tables and plots

    rod_obs_make_labels


    %% Compute observer performance metrics

    % Errors in observer state estimates
    Y_errors = repmat(Y, 1, n_obs) - Y_est;

    % Mean-squared errors
    Y_MSE = mean(Y_errors.^2, 1);

    % Calculate metrics for steady-state periods and after steps
    tau_ss = 0.43*3;  % approximate settling time
    n_settle = ceil(tau_ss/Ts);
    ss_periods = steady_state_periods(Pd, n_settle);
    Y_MSE_tr = mean(Y_errors(~ss_periods,:).^2, 1);
    Y_MSE_ss = mean(Y_errors(ss_periods,:).^2, 1);

    % Record number of samples in MSE calculations
    nT_Y_MSE = size(Y_errors, 1);
    nT_Y_MSE_tr = sum(~ss_periods);
    nT_Y_MSE_ss = sum(ss_periods);
    assert(nT_Y_MSE_tr + nT_Y_MSE_ss == nT_Y_MSE)

    % Calculate variance of estimates during steady-state periods
    trans_idxs = transition_periods(Pd);
    n_resp = numel(trans_idxs);
    sq_devs = nan(size(Y_est));
    for i = 1:n_resp
        idx = trans_idxs{i};
        idx1 = idx(1) + n_settle - 1;
        avgs = mean(Y_est(idx1:idx(2), :));
        sq_devs(idx1:idx(2), :) = (Y_est(idx1:idx(2), :) - avgs).^2;
    end
    sq_devs = sq_devs(~isnan(sq_devs(:, 1)), :);
    Y_var_ss = mean(sq_devs);
    
    % Calculate mean-squared differences
    Y_diffs = [nan(1, size(Y_est, 2)); diff(Y_est, 1)];
    Y_MSD_ss = nanmean(Y_diffs(ss_periods,:).^2, 1);

    % Display summary table
    mse_table = array2table([Y_MSE' Y_MSE_tr' Y_MSE_ss' Y_var_ss' Y_MSD_ss'], ...
        'RowNames', obs_labels, ...
        'VariableNames', metrics_labels ...
    );
    % Transpose the table (complicated in MATLAB):
    mse_table_tr = rows2vars(mse_table);
    mse_table_tr = removevars(mse_table_tr, 'OriginalVariableNames');
    mse_table_tr.Properties.RowNames = {'MSE', ...
        'MSE in transitions', 'MSE in steady-state', ...
        'Variance in steady-state', 'MSD in steady-state'};

    % Display MSE results
    mse_table_tr

    % Compute errors in MKF observer estimates (if used)
    MKF_Y_errors = cell(size(sim_out.MKF_Y_est));
    MKF_Y_MSE = cell(1, n_obs_mkf);
    for f = 1:n_obs_mkf
        obs = observers{observers_mkf(f)};
        MKF_Y_MSE{f} = size(sim_out.MKF_Y_est{f}, 1);
        % Compute errors in multiple filter state estimates
        % Note: First estimates are at k=1
        MKF_Y_errors{f} = repmat(Y, 1, obs.n_filt) - sim_out.MKF_Y_est{f};
        MKF_Y_MSE{f} = mean(MKF_Y_errors{f}(2:end, :).^2, 1);
    end


    %% Combine all parameters and results and add to summary results file

    % Model parameters
    sys.name = model_name;
    sys.A = A;
    sys.B = B;
    sys.C = C;
    sys.D = D;
    sys.Ts = Ts;
    sys_params = objects2tablerow(containers.Map({'sys'}, {sys}));

    % Observer parameters
    rv_params = objects2tablerow( ...
        containers.Map({'epsilon', 'sigma_wp', 'sigma_M'}, ...
            {epsilon, sigma_wp, sigma_M}) ...
    );
    obs_params = cell(1, n_obs);
    for f = 1:n_obs
        obs = observers{f};
        params = get_obs_params(obs);
        obs_params{f} = objects2tablerow(containers.Map({obs.label}, {params}));
    end
    obs_params = horzcat(obs_params{:});

    % Simulation settings
    sim_params = table(p_case, i_in_seq, t_stop, Ts, nT, nu, ny, n_obs);

    % Observer metrics
    obs_metrics = [table(tau_ss, nT_Y_MSE, nT_Y_MSE_tr, nT_Y_MSE_ss) ...
        array2table(reshape(mse_table.Variables', [], ...
            n_obs*n_metrics), 'VariableNames', obs_metrics_labels);
    ];

    % Summary table
    summary_results = [ ...
        array2tablerow(datetime(), 'Time') ...
        sim_params ...
        sys_params ...
        array2tablerow(obs_labels, 'obs') ...
        rv_params ...
        obs_params ...
        obs_metrics ...
    ];

    % Save to csv file
    filename = sprintf('rod_obs_sim_%d_summary.csv', p_case);
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

    % Choose length of responses in sample periods
    nT_resp = 51;  % 0 to 1.5 hrs

    % Find step changes in disturbance
    [idxs, diffs] = find_step_periods(Pd, n_settle, nT_resp);
    n_resp = numel(idxs);
    fprintf("Step responses identified: %d\n", n_resp)
    Y_resps = nan(nT_resp, n_resp);
    Y_est_resps = repmat({nan(nT_resp, n_resp*1)}, 1, n_obs);
    for i = 1:n_resp
        idx = idxs{i};
        nT_max = min(diff(idx)+1, nT_resp);
        Y_resps(:, i) = Y(idx(1):idx(1)+nT_max-1);
        for f = 1:n_obs
            Y_est_resps{f}(:, i) = Y_est(idx(1):idx(1)+nT_max-1, ...
                (f-1)*ny+1:f*ny);
        end
    end
    % Convert Y_est_resp_obs to tables
    for f = 1:n_obs
        label = strcat('yk_est_', obs_labels{f});
        Y_est_resps{f} = array2table_with_name(Y_est_resps{f}', label, '_');
    end
    resp_data = [ ...
        array2table_with_name(repmat(i_in_seq, n_resp, 1), 'i_in_seq', '_') ...
        table(diffs) ...
        array2table_with_name(Y_resps', 'yk', '_') ...
        horzcat(Y_est_resps{:}) ...
    ];

    % Save to csv file
    filename = sprintf('rod_obs_sim_%d_resps.csv', p_case);
    if isfile(fullfile(results_dir, filename))
        % Load existing results and combine
        resp_data_existing = readtable(fullfile(results_dir, filename));
        fprintf("Existing step responses loaded from file: %s\n", filename)
        resp_data = vertcat(resp_data_existing, resp_data);
    end

    % Save all results to file
    writetable(drop_empty_cols(resp_data), fullfile(results_dir, filename));
    fprintf("Step responses saved to file: %s\n", filename)

end


% For results plots, run this file next:
%rod_obs_sim_plots

fprintf("Run rod_obs_sim_plots.m to produce plots.\n")

% To calculate evaluation metrics, run this file:
%rod_obs_sim_plots

fprintf("Run rod_obs_calc_metrics.m to calculate evaluation metrics.\n")
