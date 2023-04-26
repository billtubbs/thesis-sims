% Makes adjustment to Q(2, 2) of Kalmnan filters (KF1, KF2, KF3)
%
% This script is called by run_obs_sim_spec.m
%

f = sim_spec.observers.params.f;  % fusion horizon
m = sim_spec.observers.params.m;  % maximum number of shocks
d = sim_spec.observers.params.d;  % spacing parameter

% Find MKF observer to be adjusted
i_obs = cellfun(@(x) strcmp(x.label, 'MKF_SF95'), observers);

% Multiple model filter 1
label = 'MKF_SF95';
P0 = 1000*eye(n);
switch n
    case 2
        Q0 = diag([q1 0]);
    case 4
        Q0 = diag([q1 q2 0 0]);
end
R = diag(sigma_M.^2);
observers{i_obs} = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

MKF_SF95 = observers{i_obs};
