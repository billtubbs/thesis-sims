%% Display and save simulation results
%
% Saves selected simulation results to csv files. This
% script is run by run_obs_sim_spec after each simulation.
% sim_out is a struct produced by run_simulation_obs.
%
% Before running this script make sure the following
% variables exist in the current workspace:
%  - sim_spec
%  - results_dir

% Dependencies
addpath("../data-utils/")

t = sim_out.data{:,'t'};
U = sim_out.data{:,'U'};
alpha = sim_out.data{:,'alpha'};
X = sim_out.data{:,'X'};
X_est = sim_out.data{:,'X_est'};
P = sim_out.data{:,'P'};
Y = sim_out.data{:,'Y'};
Y_m = sim_out.data{:,'Y_m'};
Y_est = sim_out.data{:,'Y_est'};
E_obs = sim_out.data{:,'E_obs'};

% Only save full simulation outputs if filename specified
if isfield(sim_spec.outputs, 'results_file')
    filename = sim_spec.outputs.results_file;
    writetable(sim_out.data, fullfile(results_dir, filename));

    % Save results from multiple model filters (if used)
    f_mkf_sp = 0;
    for f = 1:n_obs_mkf

        label = observers{observers_mkf(f)}.label;
        MKF_sim_results = [sim_out.data(:, {'k', 't'}) ...
            array2table_with_name(sim_out.MKF_i{f}, 'i', '_') ...
            array2table_with_name(sim_out.MKF_p_seq_g_Yk{f}, 'p_seq_g_Yk', '_') ...
            array2table_with_name(sim_out.MKF_X_est{f}, 'X_est', '_') ...
            array2table_with_name(sim_out.MKF_trP{f}, 'trP', '_') ...
        ];

        if startsWith(observers{observers_mkf(f)}.type, "MKF_SP")
            f_mkf_sp = f_mkf_sp + 1;
            MKF_sim_results = [MKF_sim_results ...
                array2table_with_name(sim_out.MKF_SP_f_main{f_mkf_sp}, 'f_main', '_') ...
                array2table_with_name(sim_out.MKF_SP_f_hold{f_mkf_sp}, 'f_hold', '_')
            ];
        end

        [folder, baseFileName, extension] = fileparts(filename);
        mkf_filename = sprintf("%s%s_%s%s", folder, baseFileName, label, extension);
        writetable(MKF_sim_results, fullfile(results_dir, mkf_filename));

    end

end


%% Prepare labels for tables and plots

rod_obs_sim_labels


%% Compute root-mean-squared errors and save results

% Compute errors in observer state and output estimates
X_errors = repmat(X, 1, n_obs) - X_est;
Y_errors = repmat(Y, 1, n_obs) - Y_est;

% Compute RMSEs of observer state and output estimates
X_rmse = sqrt(mean(X_errors.^2, 1));
Y_rmse = sqrt(mean(Y_errors.^2, 1));

% RMSE summary table
rmse_table = array2table([reshape(X_rmse, n, n_obs); reshape(Y_rmse, ny, n_obs)], ...
    'RowNames', [x_est_labels y_est_labels], ...
    'VariableNames', cellstr(obs_labels) ...
);
if sim_spec.setup.verbose
    % Display RMSE summary table
    fprintf("\nRMSE results\n")
    disp(rmse_table)
end

% Compute time series of MSE over simulation to check convergence
Y_crmse = cum_RMSEs(Y_errors);

% Compute errors in MKF filter estimates (if used)
MKF_X_errors = cell(size(sim_out.MKF_X_est));
MKF_X_rmse = cell(1, n_obs_mkf);
for f = 1:n_obs_mkf
    obs = observers{observers_mkf(f)};
    MKF_X_rmse{f} = size(sim_out.MKF_X_est{f}, 1);
    % Compute errors in multiple filter state estimates
    % Note: First estimates are at k=1
    MKF_X_errors{f} = repmat(X, 1, obs.nh) - sim_out.MKF_X_est{f}(:,1:obs.n*obs.nh);
    MKF_X_rmse{f} = mean(MKF_X_errors{f}.^2, 1);
end

% Compute errors immediately after shocks
n_shocks = sum(alpha, 1)';
shock_period = 10;
[rows, cols] = ind2sub(size(alpha), find(alpha == 1));
X_rmses_after_shocks = cell(1, n_obs);
for f = 1:n_obs
    sum_sq_errors = zeros(shock_period, n);
    counts = zeros(1, n);
    for i = 1:sum(n_shocks)
        errors = X_errors(rows(i):min(end, rows(i)+shock_period-1), (1:n)+(f-1)*n);
        sum_sq_errors(1:size(errors,1),:) = sum_sq_errors(1:size(errors,1),:) + errors.^2;
        counts = counts + size(errors, 1);
    end
    shock_rmses = sum_sq_errors ./ counts;
    X_rmses_after_shocks{f} = shock_rmses;
end


%% Combine all parameters and results and add to summary results file

rmse_results = flattentable(rmse_table, 'RMSE');
sys_model_params = objects2tablerow( ...
    containers.Map({'u_known', 'y_meas', 'params'}, ...
    {u_known, y_meas, params}) ...
);
input_params = objects2tablerow( ...
    containers.Map({'epsilon', 'sigma_wp', 'sigma_M'}, ...
        {epsilon, sigma_wp, sigma_M}) ...
);
obs_params = cell(1, n_obs);
for f = 1:n_obs
    obs = observers{f};
    params = get_obs_params(obs);
    params.type = obs.type;
    obs_params{f} = objects2tablerow(containers.Map({char(obs.label)}, {params}));
end
obs_params = horzcat(obs_params{:});
results_table = [ ...
    array2tablerow(datetime(), 'Time') ...
    table(sim_label, seed, t_stop, Ts, nT, nu, ny, n) ...
    sys_model_params ...
    array2tablerow(obs_labels, 'obs') ...
    input_params ...
    obs_params ...
    rmse_results ...
];

% Save to csv file
filename = sim_spec.outputs.summary_file;
if isfile(fullfile(results_dir, filename))
    % Load existing results and combine
    existing_results = readtable(fullfile(results_dir, filename), ...
        'TextType','string');
    fprintf("Existing summary file: %s\n", filename)
    results_table = outerjoin(existing_results, results_table, ...
        'MergeKeys', true);
end

% Save all results to file
writetable(results_table, fullfile(results_dir, filename));
fprintf("Summary saved to file:\n%s\n", fullfile(results_dir, filename))

return