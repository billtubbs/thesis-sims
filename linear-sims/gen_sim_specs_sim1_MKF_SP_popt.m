% Generate simulation specifications for tuning MKF_SP
% observer parameters on SISO system

clear all

addpath("../yaml")
addpath("../process-observers")

% Main folder where all simulation results will be saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim1_MKF_SP_popt";

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
nh_values = {3, 4, 5, 7, 10, 14, 19, 25};
n_min_values = {1, 2, 3, 4, 5, 6, 7, 9, 12, 16, 21};

fprintf("Generating simulation specification files...\n")

% Define struct to store info

% Basic info
data.info.label = "<placeholder>";  % this will be changed below
data.info.description = ...
    "Simulations of MKF observer with sequence pruning algorithm " + ...
    "on linear SISO system with varying covariance parameter to " + ...
    "determine best tuning.";
fprintf("%16s: %s\n", "Simulation", sim_name)

% System
data.system.setup_script = "sys_rodin_step.m";
data.system.label = "rodin_step";

% Input generation
data.inputs.setup_script = "rod_obs_sim_inputs.m";

% Observers
data.observers.setup_script = "obs_rodin_step.m";
data.observers.adj_script = "obs_rodin_step_MKF_SP_nh_nmin.m";

% Simulation setup
data.setup.plot_dir = fullfile(sims_dir, sim_name, "plots");
data.setup.results_dir = fullfile(sims_dir, sim_name, "results");
data.setup.seed = 6;
data.setup.observers = ["MKF_SP1"];
data.setup.t_stop = 2500;
data.setup.verbose = false;

% Output handling
data.outputs.setup_script = "rod_obs_sim_outputs.m";
%data.outputs.results_file = "<placeholder>";  % don't save outputs
data.outputs.summary_file = "rod_obs_sim_outputs_summary.csv";

[nh_idx, n_min_idx] = ndgrid(nh_values, n_min_values);
n_combs = numel(nh_idx);

%n_seeds = numel(seed_values);
%n_adj = numel(adj_values);
%for i = 1:n_seeds

% Load system
sys_rodin_step

n_specs = 0;
for i_comb = 1:n_combs

    % Get parameter values
    nh = nh_idx{i_comb};  % fusion horizon
    n_min = n_min_idx{i_comb};  % maximum number of shocks

    % Check combination is valid
    nw = 1;
    n_hold = nw*n_min;
    n_main = nh - n_hold;

    if n_main < nw
        fprintf("Skipping simulation due to nh - nw*n_min < nw\n")
        continue
    end
    
    n_specs = n_specs + 1;

    % Label
    data.info.label = sprintf("%s_%03d", sim_name, n_specs);
    %data.outputs.results_file = sprintf("rod_obs_sim_outputs_%02d_%02d.csv", i, j);

    % Change simulation variables
    %data.setup.seed = seed_values(i);
    data.observers.params.nh = nh;
    data.observers.params.n_min = n_min;

    % Save to Yaml file
    filename = sprintf("sim_spec_%03d.yml", n_specs);
    yaml.dumpFile(fullfile(filepath, filename), data, "block")
    fprintf("%16s: %s\n", "File created", filename)

end