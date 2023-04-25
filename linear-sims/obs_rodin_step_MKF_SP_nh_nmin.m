% Makes adjustment to Q(2, 2) of Kalmnan filters (KF1, KF2, KF3)
%
% This script is called by run_obs_sim_spec.m
%

nh = sim_spec.observers.params.nh;  % number of filters (hypotheses)
n_min = sim_spec.observers.params.n_min;  % minimum life

% Get observer label from spec file - should be 'MMKF_SP1'
assert(numel(sim_spec.setup.observers) == 1)
label = string(sim_spec.setup.observers(1));

% Find MKF observer to be adjusted
i_obs = cellfun(@(x) strcmp(x.label, label), observers);
assert(any(i_obs), "ValueError: observer not found.")

P0 = 1000*eye(n);
switch n
    case 2
        Q0 = diag([q1 0]);
    case 4
        Q0 = diag([q1 q2 0 0]);
end
R = diag(sigma_M.^2);
io.u_known = u_known;
io.y_meas = true(ny, 1);  % not currently used
observers{i_obs} = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,label);

MKF_SP1 = observers{i_obs};
