% Generate simulation specifications for tuning MKF_SP
% observer parameters on 2x2 system

clear all

addpath("~/yaml")

% Main folder where all simulation results will be saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim2_MKF_SP_popt";

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
% floor(6 + 1.5 .^ [0 1 2 3 4 5 6 7 8 9 10 11 12 13])
nh_values = {10 11 12 14 16 20 26 34 47 66 95 138 203};
% floor(1.5 .^ [0 1 2 3 4 5 6 7 8 9 10 11 12])
n_min_values = {1 2 3 5 7 11 17 25 38 57 86 129};

fprintf("Generating simulation specification files...\n")

% Define struct to store info

% Basic info
data.info.label = "<placeholder>";  % this will be changed below
data.info.description = ...
    "Simulations of MKF observer with sequence pruning algorithm " + ...
    "on linear 2x2 system with varying covariance parameter to " + ...
    "determine best tuning.";
fprintf("%16s: %s\n", "Simulation", sim_name)

% System
data.system.setup_script = "sys_rodin_step_2x2sym2.m";
data.system.label = "rodin_step_2x2sym2";

% Input generation
data.inputs.setup_script = "rod_obs_sim_inputs.m";

% Observers
data.observers.setup_script = "obs_rodin_step_2x2.m";
data.observers.adj_script = "obs_rodin_step_MKF_SP_nh_nmin.m";

% Simulation setup
data.setup.plot_dir = fullfile(sims_dir, sim_name, "plots");
data.setup.results_dir = fullfile(sims_dir, sim_name, "results");
data.setup.seed = 0;
data.setup.observers = "MKF_SP1";
data.setup.t_stop = 5000;
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

n_specs = 0;
for i_comb = 1:n_combs

    % Get parameter values
    nh = nh_idx{i_comb};  % fusion horizon
    n_min = n_min_idx{i_comb};  % maximum number of shocks

    % Check combination is valid
    nw = 2;
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