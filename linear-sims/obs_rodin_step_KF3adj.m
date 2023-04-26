% Makes adjustment to Q(2, 2) of Kalman filters (KF1, KF2, KF3)
%
% This script is called by run_obs_sim_spec.m for tuning the KF3
% Q parameters.
%

% To generate log ranges of values use:
% logspace(-2, 2, 9)
% 0.0100    0.0316    0.1000    0.3162    1.0000    3.1623   10.0000   31.6228  100.0000
% logspace(-1, 1, 9)
% 0.1000    0.1778    0.3162    0.5623    1.0000    1.7783    3.1623    5.6234   10.0000
% logspace(-1, 1, 17)
% 0.5623    0.7499    1.0000    1.3335    1.7783

% These values were used to generate the plot used in thesis report:
% 0.0100    0.0316    0.1000    0.1778    0.3162    0.5623    0.7499
% 1.0000    1.3335    1.7783    3.1623   10.0000   31.6228  100.0000

sigma_wp_adj = sim_spec.observers.params.sigma_wp;

switch KF1.n
    case 2  % SISO system model

        KF3.model.Q(2, 2) = sigma_wp_adj.^2;

    case 4  % 2x2 system model

        KF3.model.Q(3, 3) = sigma_wp_adj.^2;
        KF3.model.Q(4, 4) = sigma_wp_adj.^2;

end