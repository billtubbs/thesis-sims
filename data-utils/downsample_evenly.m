function [Xmean, Xmin, Xmax] = downsample_evenly(X, nT_min)
    if nargin < 2
        nT_min = 500;
    end
    nT = size(X, 1);
    n = size(X, 2);
    nT_new = cumprod(fliplr(factor(nT)));
    i = find((nT_new >= nT_min) & (nT_new < nT));
    if ~isempty(i)
        nT_new = nT_new(i(1));
        % X = [1 2 3 4; 5 6 7 8; 9 10 11 12]'
        % X = pagetranspose(reshape(X, 2, [], 3))
        X = reshape(X, [], nT_new, n);
        Xmean = squeeze(pagetranspose(mean(X)));
        Xmin = squeeze(pagetranspose(min(X)));
        Xmax = squeeze(pagetranspose(max(X)));
    else
        Xmean = X;
        Xmin = X;
        Xmax = X;
    end
end