% Generate simulation specifications for tuning KF
% on SISO linear system

clear all

addpath("~/yaml")

% Main folder where all simulation results will be saved
sims_dir = "simulations";

% Simulation group name
sim_name = "rod_obs_sim1_3KF_Q";

% Create main directory for these simulations
base_dir = fullfile(sims_dir, sim_name);
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Subdirectory for spec files
sim_specs_dir = "sim_specs";
filepath = fullfile(base_dir, sim_specs_dir, "queue");
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Seeds for each simulation
seed_values = 0:9;

% Default value to which adjustments are made
sigma_wp_default = 0.1;

% Use these adjustment factors to vary the parameter of interest
adj_values = [ ...
    0.0100    0.0316    0.1000    0.1778    0.3162    0.5623    0.7499 ...
    1.0000    1.3335    1.7783    3.1623   10.0000   31.6228  100.0000
];

fprintf("Generating simulation specification files...\n")

% Define struct to store info

% Basic info
data.info.label = "<placeholder>";  % this will be changed below
data.info.description = ...
    "Simulations of single Kalman filter on linear SISO system with " + ...
    "varying covariance parameter to determine best tuning.";
fprintf("%16s: %s\n", "Simulation", sim_name)

% System
data.system.setup_script = "sys_rodin_step.m";
data.system.label = "rodin_step";

% Input generation
data.inputs.setup_script = "rod_obs_sim_inputs.m";

% Observers
data.observers.setup_script = "obs_rodin_step.m";
data.observers.adj_script = "obs_rodin_step_KF3adj.m";

% Simulation setup
data.setup.plot_dir = fullfile(sims_dir, sim_name, "plots");
data.setup.results_dir = fullfile(sims_dir, sim_name, "results");
data.setup.seed = "<placeholder>";  % this will be changed below;
data.setup.no_noise = false;
data.setup.observers = ["KF1", "KF2", "KF3"];
data.setup.t_stop = 2500;
data.setup.verbose = false;

% Output handling
data.outputs.setup_script = "rod_obs_sim_outputs.m";
%data.outputs.results_file = "<placeholder>";  % don't save outputs
data.outputs.summary_file = "rod_obs_sim_outputs_summary.csv";

n_seeds = numel(seed_values);
n_adj = numel(adj_values);
for i = 1:n_seeds

    for j = 1:n_adj

        % Label
        data.info.label = sprintf("%s_%02d_%02d", sim_name, i, j);
        %data.outputs.results_file = sprintf("rod_obs_sim_outputs_%02d_%02d.csv", i, j);

        % Change simulation variables
        data.setup.seed = seed_values(i);
        data.observers.params.sigma_wp = sigma_wp_default .* adj_values(j);

        % Save to Yaml file
        filename = sprintf("sim_spec_%02d_%02d.yml", i, j);
        yaml.dumpFile(fullfile(filepath, filename), data, "block")
        fprintf("%16s: %s\n", "File created", filename)

    end
end