% Makes adjustment to Q(2, 2) of Kalmnan filters (KF1, KF2, KF3)
%
% This script is called by run_obs_sim_spec.m
%

nf = sim_spec.observers.params.nf;  % fusion horizon
m = sim_spec.observers.params.m;  % maximum number of shocks
d = sim_spec.observers.params.d;  % spacing parameter

% Find MKF observer to be adjusted
i_obs = cellfun(@(x) strcmp(x.label, 'MKF_SF1'), observers);

% Multiple model filter 1
label = 'MKF_SF1';
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
observers{i_obs} = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
        Q0,R,nf,m,d,label);

MKF_SF1 = observers{i_obs};
