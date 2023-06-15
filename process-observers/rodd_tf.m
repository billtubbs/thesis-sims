function Hd = rodd_tf(theta0, Phi, d, Ts)
% rodd_tf returns the discrete transfer function of a
% randomly-occurring deterministic disturbance process as
% defined in section 2.1 of MacGregor et al. (1984).
%
% Arguments
%  theta0 : direct transmission gain
%  Phi : polynomial factor to include in denominator
%  d : number of integrators to include in the denominator
%  Ts : sample time
%
% Examples
% >> Hd = rodd_tf(1, 1, 1, 1)  % step process
% 
% Hd =
%  
%     1
%   -----
%   z - 1
%  
% Sample time: 1 seconds
% Discrete-time transfer function.
% 
% >> Hd = rodd_tf(1, 1, 2, 1)  % ramp process
% 
% Hd =
%  
%         1
%   -------------
%   z^2 - 2 z + 1
%  
% Sample time: 1 seconds
% Discrete-time transfer function.
% 
% >> Hd = rodd_tf(1, [1 -0.7], 1, 1)  % exponential process
% 
% Hd =
%  
%           1
%   -----------------
%   z^2 - 1.7 z + 0.7
%  
% Sample time: 1 seconds
% Discrete-time transfer function.
% 
    for i=1:d
        Phi = conv(Phi, [1 -1]);
    end
    Hd = tf(theta0, Phi, Ts);
end