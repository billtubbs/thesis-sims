function sys = ss_concat(s1, s2)
% sys = ss_concat(s1, s2)
% Concatenate two state-space systems. This should have
% the same effect as sys = series(s1, s2) although the 
% system matrices may not be identical
%
% Example
% >> s1 = ss(tf(1,[0.5 1],1));
% >> s2 = ss(tf([1 0],[1 -1],1));
% >> ss_concat(s1, s2)
% 
% ans =
%  
%   A = 
%        x1  x2
%    x1  -2   0
%    x2   1   1
%  
%   B = 
%        u1
%    x1   2
%    x2   0
%  
%   C = 
%        x1  x2
%    y1   1   1
%  
%   D = 
%        u1
%    y1   0
%  
% Sample time: 1 seconds
% Discrete-time state-space model.
%
    [n1, nu1, ny1, Ts1, d1] = check_model(s1);
    if ~d1
        s1.D = zeros(ny1, nu1);
    end
    [n2, nu2, ny2, Ts2, d2] = check_model(s2);
    if ~d2
        s2.D = zeros(ny2, nu2);
    end
    assert(Ts1 == Ts2, "ValueError: systems have different sample periods")
    A = [s1.A zeros(n1, n2); s2.B*s1.C s2.A];
    B = [s1.B; s2.B*s1.D];
    C = [s2.D*s1.C s2.C];
    D = s2.D*s1.D;
    sys = ss(A, B, C, D, Ts1);
end