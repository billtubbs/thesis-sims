function n_filt = n_filters(m, f, nw)
% n_filt = n_filters(m, n, nw) returns the number of
% filters needed to track combinations of up to m
% randomly-occurring shocks in f periods where nw is
% the dimension of the disturbance vector w(k).
    if nargin == 2
        nw = 1;
    end
    n_filt = 0;
    for i = 0:m
        n_filt = n_filt + nchoosek(f*nw, i);
    end
end