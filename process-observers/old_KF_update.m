% Update steps described in Zhao et al.

% Model prediction of state estimates in current 
% timestep based on previous estimates and current
% input (priori state estimate)
obj.xk_est = obj.A * obj.xk_est + obj.B * uk;

% Error covariance of state prediction
obj.P = obj.A * obj.P * A' + obj.Q;

% Error covariance of output prediction
S = obj.C * obj.P * obj.C' + obj.R;

% Update correction gain
obj.Kf = obj.P * obj.C' / S;

% Update of state estimates using measurements from
% current time step (a posteriori state estimate)
obj.xk_est = obj.xk_est + obj.Kf * (yk - obj.C * obj.xk_est);

% Update error covariance of state estimates
obj.P = obj.P - obj.Kf * obj.S * obj.Kf';

% Output estimate
obj.yk_est = obj.C * obj.xk_est;