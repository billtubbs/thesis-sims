% Script to run a set of observer simulations defined in
% Yaml text files in a specified folder.
%

clear all
addpath('~/process-observers/')

% Main folder where all simulation results are saved
sims_dir = 'simulations';

% Select which set of simulations to run

% SISO linear system (see sys_rodin_step.m)
%sim_name = "rod_obs_sim1_3KF_Q";  % tuning of Kalman filters done 2022-11-21
%sim_name = "rod_obs_sim1_3KF_seed";
%sim_name = "rod_obs_sim1_MKF_SF_popt";  % re-run 2022-12-08 after fixing MKF_SF_DI
%sim_name = "rod_obs_sim1_MKF_SF95_popt";  % re-done 2022-11-29
%sim_name = "rod_obs_sim1_MKF_SP_popt";  % re-done 2022-11-29
%sim_name = "rod_obs_sim1_all_seed";

% 2x2 linear system (see rodin_step_2x2sym2.m)
%sim_name = "rod_obs_sim2_3KF_Q";  % tuning of Kalman filters
%sim_name = "rod_obs_sim2_MKF_SF_popt";  % re-ran 2022-12-08 after fixing MKF_SF_DI
sim_name = "rod_obs_sim2_MKF_SF95_popt";  % re-done 2023-04-24
%sim_name = "rod_obs_sim2_MKF_SP_popt";  % re-done 2023-04-24 with new param values
%sim_name = "rod_obs_sim2_3KF_seed";
%sim_name = "rod_obs_sim2_all_seed";


% Subdirectory for spec files
sim_specs_dir = 'sim_specs';

filepath = fullfile(sims_dir, sim_name, sim_specs_dir, 'queue');
file_pattern = "*.yml";
files_info = dir(fullfile(filepath, file_pattern));

if ~isfolder(fullfile(sims_dir, sim_name, sim_specs_dir, 'done'))
    mkdir(fullfile(sims_dir, sim_name, sim_specs_dir, 'done'))
end

% Process all files in those folders.
n_files = length(files_info);
if n_files == 0
    fprintf("No sim spec files found in queue folder:\n%s\n", filepath)
else
    fprintf("\nStarting multiple simulations...\n")
    fprintf("%16s: %d\n", "No. of files found", n_files)
end

start_time = datetime();
fprintf("%16s: %s\n", "Start time", string(start_time,"HH:mm:ss"))
for i = 1:n_files
    filepath = files_info(i).folder;
    filename = files_info(i).name;
    fprintf("\nSimulation %d of %d\n", i, n_files)

    % run simulation catching any errors
    try
        run_obs_sim_spec(filename, filepath)
    catch ME
        fprintf("Simulation failed due to the following error:\n")
        disp(getReport(ME))
        move_from = fullfile(filepath, filename);
        move_to = fullfile(sims_dir, sim_name, sim_specs_dir, 'failed');
        if ~isfolder(move_to)
            mkdir(move_to);
        end
        movefile(move_from, move_to)
        continue
    end

    move_from = fullfile(filepath, filename);
    move_to = fullfile(sims_dir, sim_name, sim_specs_dir, 'done');
    movefile(move_from, move_to)

end
end_time = datetime();
fprintf("%16s: %s\n", "End time", string(end_time,"HH:mm:ss"))
fprintf("%16s: %s\n", "Duration", string(end_time - start_time))
