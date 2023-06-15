% Generate simulation specifications for simulations
% on SISO linear system

clear all

addpath("../yaml")

% Main folder where all simulation results will be saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim1_all_seed";

% Create main directory for these simulations
base_dir = fullfile(sims_dir, sim_name);
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Subdirectory for spec files
sim_specs_dir = 'sim_specs';
filepath = fullfile(base_dir, sim_specs_dir, 'queue');
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Seeds for each simulation
seed_values = 0:10;

% Select simulation(s) for which to run a second simulation
% without measurement noise
seed_values_no_noise = 6;

fprintf("Generating simulation specification files...\n")

% Define struct to store info

% Basic info
data.info.label = "<placeholder>";  % this will be changed below
data.info.description = ...
    "Simulations of all observers on SISO system with " + ...
    "different random seeds.";
fprintf("%16s: %s\n", "Simulation", sim_name)

% System
data.system.setup_script = "sys_rodin_step.m";
data.system.label = "rodin_step";

% Input generation
data.inputs.setup_script = "rod_obs_sim_inputs.m";

% Observers
data.observers.setup_script = "obs_rodin_step_opt.m";

% Simulation setup
data.setup.plot_dir = fullfile(sims_dir, sim_name, "plots");
data.setup.results_dir = fullfile(sims_dir, sim_name, "results");
data.setup.seed = "<placeholder>";  % this will be changed below;
data.setup.observers = ["KF1", "KF2", "KF3", "SKF", "MKF_SF95", ...
    "MKF_SF95_2", "MKF_SF1", "MKF_SP1", "MKF_SP2"];
data.setup.t_stop = 2500;
data.setup.verbose = false;

% Output handling
data.outputs.setup_script = "rod_obs_sim_outputs.m";
data.outputs.results_file = "<placeholder>";  % this will be changed below;
data.outputs.summary_file = "rod_obs_sim_outputs_summary.csv";

n_seeds = numel(seed_values);
for i = 1:n_seeds

    % Label
    data.info.label = sprintf("%s_%03d", sim_name, i);
    data.outputs.results_file = sprintf("rod_obs_sim_outputs_%03d.csv", i);

    % Change variables
    data.setup.seed = seed_values(i);
    if ismember(seed_values(i), seed_values_no_noise)
        data.setup.no_noise = true;
    else
        data.setup.no_noise = false;
    end

    % Save to Yaml file
    filename = sprintf("sim_spec_%03d.yml", i);
    yaml.dumpFile(fullfile(filepath, filename), data, "block")
    fprintf("%16s: %s\n", "File created", filename)

end
