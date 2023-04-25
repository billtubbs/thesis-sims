% Generate simulation specifications for tuning MKF_SF95
% observer parameters on second 2x2 linear system

clear all

addpath("~/yaml")
addpath("~/process-observers/")

% Main folder where all simulation results will be saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim2_MKF_SF95_popt";

% Create main directory for these simulations
base_dir = fullfile(sims_dir, sim_name);
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Subdirectory for spec files
sim_specs_dir = 'sim_specs';
filepath = fullfile(base_dir, sim_specs_dir, "queue");
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Parameter values to search
n_di_values = {3, 5, 10, 15, 20};  % number of detection intervals
m_values = {1, 2, 3};
d_values = {1, 2, 3, 5, 10};

fprintf("Generating simulation specification files...\n")

% Define struct to store info

% Basic info
data.info.label = "<placeholder>";  % this will be changed below
data.info.description = ...
    "Simulations of MKF observer with sequence fusion algorithm " + ...
    "(1995) on linear 2x2 system with varying covariance parameter to " + ...
    "determine best tuning.";
fprintf("%16s: %s\n", "Simulation", sim_name)

% System
data.system.setup_script = "sys_rodin_step_2x2sym2.m";
data.system.label = "rodin_step_2x2sym2";

% Input generation
data.inputs.setup_script = "rod_obs_sim_inputs.m";

% Observers
data.observers.setup_script = "obs_rodin_step_2x2.m";
data.observers.adj_script = "obs_rodin_step_MKF_SF95_fmd.m";

% Simulation setup
data.setup.plot_dir = fullfile(sims_dir, sim_name, "plots");
data.setup.results_dir = fullfile(sims_dir, sim_name, "results");
data.setup.seed = 0;
data.setup.observers = ["MKF_SF95"];
data.setup.t_stop = 5000;
data.setup.verbose = false;

% Output handling
data.outputs.setup_script = "rod_obs_sim_outputs.m";
%data.outputs.results_file = "<placeholder>";  % don't save outputs
data.outputs.summary_file = "rod_obs_sim_outputs_summary.csv";

[n_di_idx, m_idx, d_idx] = ndgrid(n_di_values, m_values, d_values);
n_combs = numel(n_di_idx);

nh_max = 200;
beta_min = 0.85;

% Load system
sys_rodin_step_2x2sym2

n_specs = 0;
for i_comb = 1:n_combs

    % Get parameter values
    n_di = n_di_idx{i_comb};  % fusion horizon
    m = m_idx{i_comb};  % maximum number of shocks
    d = d_idx{i_comb};  % spacing parameter
    f = n_di * d;  % fusion horizon

    % Construct multiple model observer - to check beta and nh_max
    q1 = 0.01; q2 = 0.01;
    label = 'MKF_SF95';
    P0 = 1000*eye(n);
    Q0 = diag([q1 q2 0 0]);
    R = diag(sigma_M.^2);
    io.u_known = u_known;
    io.y_meas = y_meas;
    MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
        Q0,R,f,m,d,label);
    assert(MKF_SF95.f == f)

    if MKF_SF95.nh_max > nh_max
        fprintf("Skipping simulation due to nh_max > %d\n", nh_max)
        continue
    elseif MKF_SF95.beta < beta_min
        fprintf("Skipping simulation due to beta < %.2f\n", beta_min)
        continue
    end

    n_specs = n_specs + 1;

    % Label
    data.info.label = sprintf("%s_%03d", sim_name, n_specs);
    %data.outputs.results_file = sprintf("rod_obs_sim_outputs_%02d_%02d.csv", i, j);

    % Change simulation variables
    %data.setup.seed = seed_values(i);
    data.observers.params.f = f;
    data.observers.params.m = m;
    data.observers.params.d = d;
    
    % Save to Yaml file
    filename = sprintf("sim_spec_%03d.yml", n_specs);
    yaml.dumpFile(fullfile(filepath, filename), data, "block")
    fprintf("%16s: %s\n", "File created", filename)

end
