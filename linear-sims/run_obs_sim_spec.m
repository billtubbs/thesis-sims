function run_obs_sim_spec(filename, filepath)
% run_obs_sim_spec(filename, filepath)
% Sets up and runs one simulation experiment based on
% the settings specified in the Yaml file with the
% specified filepath. If the filepath argument is not 
% provided, it looks in the current directory.
%

    if nargin < 2
        filepath = '';
    end

    addpath('~/yaml')

    % Load simulation specification file
    fprintf("%16s: '%s'\n", "Spec. file", filename)
    sim_spec = yaml.loadFile(fullfile(filepath, filename));
    sim_label = sim_spec.info.label;
    fprintf("%16s: %s\n", "Label", sim_label)

    % Prepare sub-directories if they don't already exist
    plot_dir = sim_spec.setup.plot_dir;
    results_dir = sim_spec.setup.results_dir;
    if ~isfolder(plot_dir)
        fprintf("%16s: %s\n", "New sub-directory", plot_dir)
        mkdir(plot_dir);
    end
    if ~isfolder(results_dir)
        fprintf("%16s: %s\n", "New sub-directory", results_dir)
        mkdir(results_dir);
    end

    % Set random number generator seed
    seed = sim_spec.setup.seed;
    rng(seed)

    % Load system from file
    run(sim_spec.system.setup_script)
    sys_label = sim_spec.system.label;
    fprintf("%16s: %s\n", "System", sys_label)

    % Simulation setup
    t_stop = sim_spec.setup.t_stop;
    nT = ceil(t_stop / Ts);
    fprintf("%16s: %d\n", "Time steps", nT)

    % Time vector
    t = Ts * (0:nT)';

    % Load observers from file
    run(sim_spec.observers.setup_script)
    % NOTE: this file should produce a cell array of
    % observer objects called 'observers':
    assert(numel(observers) > 0)

    % Run script to make adjustments to system parameters
    if ismember('params', fieldnames(sim_spec.system))
        if ismember('adj_script', fieldnames(sim_spec.system))
            fprintf("Adjusting system parameters...\n")
            run(sim_spec.system.adj_script)
        end
    end

    % Generate input data
    run(sim_spec.inputs.setup_script)  % this also adds the SKF observer
    fprintf("%16s: %s\n", "Inputs", strjoin(...
        input_data.Properties.VariableNames, ', '))
    assert(size(U, 1) == nT+1)
    assert(size(V, 1) == nT+1)
    assert(size(W, 1) == nT+1)

    all_obs_labels = cellfun(@(x) char(x.label), observers, 'UniformOutput', false);
    assert(numel(unique(all_obs_labels)) == numel(all_obs_labels))

    % Run script to make adjustments to observers
    if ismember('params', fieldnames(sim_spec.observers))
        if ismember('adj_script', fieldnames(sim_spec.observers))
            fprintf("Adjusting observer parameters...\n")
            run(sim_spec.observers.adj_script)
        end
    end

    % Observers to include in simulation
    obs_labels = string(sim_spec.setup.observers);
    n_obs = numel(obs_labels);
    fprintf("%16s: %s\n", "Observers", strjoin(obs_labels, " "))

    % Prepare cell array of selected observers
    obs_sel = cell(1, n_obs);
    if ~isfield(sim_spec.setup, "no_noise")
        sim_spec.setup.no_noise = false;
    end
    if sim_spec.setup.no_noise
        obs_sel_no_noise = cell(1, n_obs);
    end
    for i = 1:n_obs
        matches = strcmp(all_obs_labels, obs_labels{i});
        assert(sum(matches) > 0, "Observer label not found.")
        assert(sum(matches) == 1, "Duplicate observer labels.")
        if any(matches)
            obs_sel{i} = observers{matches};
            if sim_spec.setup.no_noise
                obs_sel_no_noise{i} = copy(observers{matches});
            end
        else
            error('Observer "%s" not found.', obs_labels{i})
        end
    end
    % Keep only selected obervers
    observers = obs_sel;

    % Identify multi-model observers
    n_obs_mkf = 0;
    observers_mkf = double.empty(1, 0);
    for i = 1:n_obs
        if startsWith(observers{i}.type, "MKF")
            n_obs_mkf = n_obs_mkf + 1;
            observers_mkf(n_obs_mkf) = i;
        end
    end
    fprintf("%16s: %d\n", "No. of MKFs", n_obs_mkf)

    % Prepare transition functions for simulating system
    [state_fcn, meas_fcn, params] = make_trans_funcs_ssd(A, B, C, D);

    % Simulate whole system and observers
    fprintf("Starting simulation...\n")
    sim_out = run_simulation_obs(Ts, input_data, state_fcn, ...
        meas_fcn, params, x0, u_known, observers);
    fprintf("Simulation complete.\n")
    fprintf("%16s: (%d, %d)\n", "Outputs size", size(sim_out.data))
    fprintf("%16s: %s\n", "Outputs", strjoin(...
        sim_out.data.Properties.VariableNames, ', '))
    if sim_spec.setup.verbose
        % Display simulation outputs
        sim_out.data
    end

    % Run script to save outputs
    results_dir = sim_spec.setup.results_dir;
    run(sim_spec.outputs.setup_script)

    if sim_spec.setup.no_noise
        % Run no-noise simulation if selected
        sim_out_with_noise = sim_out;
        [state_fcn, meas_fcn, params] = make_trans_funcs_ssd(A, B, C, D);
        fprintf("Starting 'no noise' simulation...\n")
        input_data{:, 'V'} = 0;  % set measurement noise to zero
        sim_out = run_simulation_obs(Ts, input_data, state_fcn, ...
            meas_fcn, params, x0, u_known, obs_sel_no_noise);
        fprintf("Simulation complete.\n")
        fprintf("%16s: (%d, %d)\n", "Outputs size", size(sim_out.data))
        fprintf("%16s: %s\n", "Outputs", strjoin(...
            sim_out.data.Properties.VariableNames, ', '))    
        if sim_spec.setup.verbose
            sim_out.data
        end
        % Run script to save no-noise outputs
        results_dir = sim_spec.setup.results_dir + "_no_noise";
        if ~isfolder(results_dir)
            mkdir(results_dir)
        end
        run(sim_spec.outputs.setup_script)
        % Restore both sets of simulation results
        sim_out_no_noise = sim_out;
        sim_out = sim_out_with_noise;
    end

    % Run script to make plots (if selected)
    if isfield(sim_spec, "plots")
        plot_dir = sim_spec.setup.plot_dir;
        run(sim_spec.plots.setup_script)
    end

end