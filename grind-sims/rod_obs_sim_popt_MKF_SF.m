% Simulate the MKF_SF_RODD95 observer on simulated measurement 
% data from grinding simulation model with a range of different
% parameter settings for optimization.
%
% Input files:
%  - rod_obs_P2DcTd4.m - process model and observers
%

clear all

% Specify path to observer functions
addpath("../process-observers")
addpath("../data-utils")
addpath("../plot-utils")

% Sub-directories used
data_dir = 'data';
results_dir = 'results';
plot_dir = 'plots';
if ~isfolder(results_dir)
    mkdir(results_dir);
end
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Specify which simulation case
p_case = 1;  % Not currently used

% Specify which data set(s) to run simulations with
% I used:
% - 1 for process model estimation (Fig. 4 in paper)
% - 2 for process model validation (model selection)
% - 3 for initial observer test (Fig. 5 in paper)
% - 6 for observer parameter optimization
% - 6 to 15 for observer Monte Carlo simulations.
i_in_seq = 7;

% Labels to identify results file
%obs_label = "MKF_SF95";
obs_label = "MKF_SF1";
sim_label = "popt_" + obs_label;

% Load observers
%rod_obs_P1Dcd4
%rod_obs_P1
%rod_obs_P1DcD5
%rod_obs_P2U
rod_obs_P2DcTd4  % observers used in IFAC paper
%rod_obs_oe125

% Generate the simulation data with the following script which 
% runs the Simulink model:
%  - sim_experiment_ore_switching.m

% Define parameter ranges
% TODO: Add 7 because it could be good
nf_values = {1, 2, 3, 4, 5, 6, 8, 10};
m_values = {1, 2, 3};
d_values = {3, 4, 5, 6, 8, 10, 12, 15};

% Determine all possible combinations
[ndi_idx, m_idx, d_idx] = ndgrid(nf_values, m_values, d_values);

% Eliminate some combinations
fprintf("Selecting parameter value combinations...\n")
n_combs = numel(ndi_idx);
selection = true(numel(ndi_idx), 1);
for i_comb = 1:n_combs
    % Parameter values
    nf = ndi_idx{i_comb};  % number of detection intervals
    m = m_idx{i_comb};  % maximum no. of shocks in fusion horizon
    d = d_idx{i_comb};  % length of detection interval

    % Some combinations are invalid
    if m > nf
        selection(i_comb) = false;
    end

    % Long fusion horizon
    f = nf * d;  % fusion horizon
    if f > 100
        selection(i_comb) = false;
    end
    %fprintf("(%d, %d, %d)\n", f, m, d)
end

% Remove unwanted combinations
ndi_idx = ndi_idx(selection);
m_idx = m_idx(selection);
d_idx = d_idx(selection);
fprintf("%d of %d possible combinations will be tested\n", ...
    numel(ndi_idx), n_combs)
n_combs = numel(ndi_idx);

% Start simulations
for i_comb = 1:n_combs

    % Create observer with parameter values
    nf = ndi_idx{i_comb};  % number of detection intervals
    m = m_idx{i_comb};  % maximum no. of shocks in fusion horizon
    d = d_idx{i_comb};  % length of detection interval
    f = nf * d;  % fusion horizon

    % Choose the observer to simulate
    i_obs = find(cellfun(@(obs) strcmp(obs.label, obs_label), observers));
    assert(numel(i_obs) == 1)
    obs = observers{i_obs};
    assert(ismember(obs.type, ["MKF_SF_RODD95" "MKF_SF_RODD"]))

    % Re-initialize observer - sequence fusion
    switch obs.type
        case "MKF_SF_RODD95"
            obs = MKFObserverSF_RODD95(model,io,obs.P0,obs.epsilon, ...
                obs.sigma_wp,obs.Q0,obs.R,f,m,d,obs.label);
        case "MKF_SF_RODD"
            obs = MKFObserverSF_RODD(model,io,obs.P0,obs.epsilon, ...
                obs.sigma_wp,obs.Q0,obs.R,nf,m,d,obs.label);
    end

    % Reject combination if it does not meet criteria
    nh_max = 200;
    beta_min = 0.85;
    if (obs.nh_max > nh_max) || (obs.beta < beta_min)
        fprintf("Rejecting due to nh_max > %d (%d)\n", nh_max, obs.nh_max)
        continue
    end

    observers = {obs};

    fprintf("\nObserver simulation %d of %d with \n", i_comb, n_combs)
    fprintf("f: %d, m: %d, d: %d, nh: %d, Input seq.: #%d\n", f, m, d, ...
        obs.nh_max, i_in_seq)

    % Load system simulation results
    if i_in_seq < 6
        nT = 300;
    else
        nT = 2460;
    end
    filename = sprintf('sim_OL_rc_est_mix_factor_%d_%d_ident.csv', nT, i_in_seq);
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
    % disturbance
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

    % Run simulation
    input_data = table(U, alpha, gamma, Pd, Y, Y_m);
    sim_out = run_obs_simulation(Ts, input_data, observers);
    observers = sim_out.observers;  % Updated observers

    
    %% Display and save simulation results
    
    % Remove semi-colon to display results table
    sim_out.data;

    % No real need to save the results
    %filename = sprintf('rod_obs_sim_%s_%d_%03d.csv', sim_label, p_case, i_comb);
    %writetable(sim_out.data, fullfile(results_dir, filename));
    %fprintf("Observer simulation results saved to file: %s\n", filename)

    % Count number of observers and MKF/AFMM observers
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
        writetable(MKF_sim_results, fullfile(results_dir, filename));
        fprintf("MKF simulation results saved to file: %s\n", filename)

    end


    %% Prepare labels for tables and plots

    rod_obs_make_labels


    %% Compute observer performance metrics

    % Approximate settling time (was 0.43*3)
    tau_ss = 1.2;

    [metrics, metrics_params, errors, metrics_labels] = ...
        calculate_obs_metrics(Y, Y_est, obs_labels, Pd, Ts, tau_ss);

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
    filename = sprintf('rod_obs_sim_%s_%d_summary.csv', sim_label, p_case);
    if isfile(fullfile(results_dir, filename))
        % Load existing results and combine
        summary_results_existing = readtable(fullfile(results_dir, filename));
        fprintf("Existing results loaded from file: %s\n", filename)
        summary_results = vertcat(summary_results_existing, summary_results);
    end

    % Save all results to file
    writetable(summary_results, fullfile(results_dir, filename));
    fprintf("Summary results saved to file: %s\n", filename)

end


% To plot results of popt run this script

%rod_obs_sim_popt_MKF_SF_plots

fprintf("run rod_obs_sim_popt_MKF_SF_plots.m to produce plots.\n")



