function [X, Y, Ym] = run_simulation_sys(models,U,V,Gamma,nT,x0)
% [X, Y, Ym] = run_simulation_sys(models,U,V,Gamma,nT,x0)
% Simulate switching system - used by the following test 
% scripts:
%
%  - test_MKFObservers_JS.m
%  - test_MKFObservers_2x2.m
%

    arguments
        models (1, :) cell
        U double
        V double
        Gamma double {mustBeInteger}
        nT (1, 1) double {mustBeInteger}
        x0 (:, 1) double = [];
    end

    [nj, n, nu, ny, Ts] = check_models(models);

    if isempty(x0)
        x0 = zeros(n, 1);
    end

    X = zeros(nT+1,n);
    Y = zeros(nT+1,ny);
    Ym = zeros(nT+1,ny);
    xk = x0;

    for i = 1:nT+1

        % Switch system
        j = Gamma(i) + 1;

        % Inputs
        uk = U(i,:)';

        % Compute y(k)
        yk = models{j}.C * xk;
        yk_m = yk + V(i);

        % Store results
        X(i, :) = xk';
        Y(i, :) = yk';
        Ym(i, :) = yk_m';

        % Compute x(k+1)
        xk = models{j}.A * xk + models{j}.B * uk;

    end
end