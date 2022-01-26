% Prepare identification data from simulation results for
% use in System Identification toolbox.

% Generate the simulation data using the following:
%  - sim_experiment_ore_switching.m

% Directory where simulation results are stored
results_dir = 'data';

% Specify which data set to load (also determines random 
% seed used for generating measurement noise)
i_sim = 15;
sim_files = { ...
    ... 1 for estimation data (estimation of models)
    'sim_OL_rc_est_mix_factor_1200_01.csv', ...
    ... 2 for validation data (model selection)
    'sim_OL_rc_est_mix_factor_1200_02.csv', ...
    ... 3 for testing
    'sim_OL_rc_est_mix_factor_1200_03.csv', ...
    ... 4-8 longer sequences for Monte-Carlo tests
    'sim_OL_rc_est_mix_factor_1200_04.csv', ...
    'sim_OL_rc_est_mix_factor_1200_05.csv', ...
    'sim_OL_rc_est_mix_factor_7680_06.csv', ...
    'sim_OL_rc_est_mix_factor_7680_07.csv', ...
    'sim_OL_rc_est_mix_factor_7680_08.csv', ...
    'sim_OL_rc_est_mix_factor_7680_09.csv', ...
    'sim_OL_rc_est_mix_factor_7680_10.csv' ...
    'sim_OL_rc_est_mix_factor_7680_11.csv', ...
    'sim_OL_rc_est_mix_factor_7680_12.csv', ...
    'sim_OL_rc_est_mix_factor_7680_13.csv', ...
    'sim_OL_rc_est_mix_factor_7680_14.csv', ...
    'sim_OL_rc_est_mix_factor_7680_15.csv' ...
};
filename = sim_files{i_sim};

% Different random number generator seed in each case
rng(i_sim, 'twister')

% Load simulation results
sim_data = readtable(fullfile(results_dir, filename));
t_data =sim_data.t;
% Check sampling period of data
Ts_data = 1/60;
assert(abs(diff(t_data(1:2)) - Ts_data) < 1e-6)
assert((size(t_data, 1) - 1) * Ts_data - t_data(end) == 0)

% Re-sample data for observers
Ts = 3/60;  % desired sample period
d = Ts / Ts_data;
assert(mod(d, 1) == 0);  % must be a whole number of samples

% Define input and output variables for model
mvs = {'BASE_ORE_MIX'};
cvs = {'SAG_J', 'SAG_POW', 'SAG_OF_P80'};
nu = numel(mvs);
ny = numel(cvs);

% Select data for identification
selected_data = sim_data(1:d:end, [{'t'} mvs cvs]);
t = selected_data.t;

% Strip initial period before equilibrium reached
t_start = 5;
selected_data = selected_data(t >= t_start, :);
nT = size(selected_data, 1) - 1;
selected_data.t = selected_data.t - t_start;
t = selected_data.t;

% Load normal operating points of process variables
op_data

% Add measurement noise to outputs
sigma_M = [1, 2, 5];
cvs_m = {'SAG_J_M', 'SAG_POW_M', 'SAG_OF_P80_M'};

V = randn(nT+1, ny) .* sigma_M;
selected_data{:, cvs_m} = selected_data{:, cvs} + V;

figure(1); clf
make_io_plots(selected_data, mvs, cvs, cvs_m)


%% Remove operating points

% Estimate initial conditions
i_steps = find(diff(selected_data{:, mvs{1}}));
y0 = mean(selected_data{1:i_steps(1), cvs});
y0_m = mean(selected_data{1:i_steps(1), cvs_m});

% Remove operating points (initial conditions)
assert(nu == 1 & isequal(mvs, {'BASE_ORE_MIX'}))
ident_data = selected_data;  % make a copy
ident_data{:, mvs} = selected_data{:, mvs} - op_pts.('mix_factor');
ident_data{:, cvs} = selected_data{:, cvs} - y0;
ident_data{:, cvs_m} = selected_data{:, cvs_m} - y0_m;

figure(2); clf
make_io_plots(ident_data, mvs, cvs, cvs_m)

u = ident_data{:, mvs};
y = ident_data{:, cvs};
y_m = ident_data{:, cvs_m};


%% Save identification dataset
filename = sprintf('sim_OL_rc_est_mix_factor_%d_%02d_ident.csv', nT, i_sim);
writetable(ident_data, fullfile(results_dir, filename))
fprintf("Identification data saved to file '%s'.\n", filename)


%% Make IDDATA objects for use in System Identification toolbox

% Identification data - filling and power only
id_data_12 = iddata(y(:, [1 2]), u, Ts);
id_data_12.TimeUnit = 'minutes';
id_data_12_m = iddata(y_m(:, [1 2]), u, Ts);
id_data_12_m.TimeUnit = 'minutes';

% Identification data - P80 only
id_data_3 = iddata(y(:, 3), u, Ts);
id_data_3.TimeUnit = 'minutes';
id_data_3_m = iddata(y_m(:, 3), u, Ts);
id_data_3_m.TimeUnit = 'minutes';

fprintf("id_data objects for sim #%d created.\n", i_sim)

% Now run system identification tool
%systemIdentification


function make_io_plots(data, mvs, cvs, cvs_m)
    nu = numel(mvs);
    ny = numel(cvs);
    t = data.t;
    n_plots = ny + nu;
    axs = nan(1, n_plots);

    for i = 1:numel(cvs)
        axs(i) = subplot(n_plots, 1, i);
        plot(t, data{:, cvs_m(i)}, '.'); hold on
        plot(t, data{:, cvs(i)}, 'Linewidth', 2)
        ylabel(cvs(i), 'Interpreter', 'none');
        grid on
    end

    for i = 1:numel(mvs)
        axs(numel(cvs) + i) = subplot(n_plots, 1, numel(cvs) + i);
        Y = data{:, mvs(i)};
        plot(t, Y, 'Linewidth', 2)
        min_max = [min(Y)-0.02 max(Y)+0.02];
        ylim(min_max);
        ylabel(mvs(i), 'Interpreter', 'none');
        grid on
        if i == numel(mvs)
            xlabel('Time')
        end
    end
    linkaxes(axs, 'x');
end

