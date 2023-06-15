function [rkm1, rk] = mode_transitions_all(nj)
% [rkm1, rk] = mode_transitions_all(nj)
% Returns two vectors representing all possible mode
% transitions between time k - 1 and k for a system
% with nj modes.
%
% Example:
% >> [rkm1, rk] = mode_transitions_all(2)
% 
% rkm1 =
% 
%      1
%      1
%      2
%      2
% 
% 
% rk =
% 
%      1
%      2
%      1
%      2
%

    [rkm1, rk] = meshgrid(1:nj, 1:nj);
    rkm1 = reshape(rkm1, [], 1);
    rk = reshape(rk, [], 1);

end