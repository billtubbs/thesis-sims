% Test model parameter-checking functions
% 
%  - check_dimensions.m
%  - check_model.m
%  - check_models.m
%  - check_model_NL.m
%  - check_models_NL.m
%

clear all
rng(0)

% Dimensions of test systems
p.n = 4; p.ny = 2; p.nu = 3;

% Sampling period
p.Ts = 0.5;

% Number of systems
nj = 3;


%% Check linear models

% Create random linear systems with given dimensions
n_models = nj;
models = cell(1, n_models);
for i = 1:n_models
    models{i} = drss(p.n,p.ny,p.nu);
    models{i}.Ts = p.Ts;
end

model = models{1};
[test.n, test.nu, test.ny] = check_dimensions(model.A, model.B, model.C);
assert(isequal([test.n, test.nu, test.ny], [p.n, p.nu, p.ny]))

[test.n, test.nu, test.ny] = check_dimensions(model.A, model.B, model.C, ...
    model.D);
assert(isequal([test.n, test.nu, test.ny], [p.n, p.nu, p.ny]))

[test.n, test.nu, test.ny, test.Ts, direct] = check_model(model);
assert(isequal(test, p))

[test_nj, test.n, test.nu, test.ny, test.Ts] = check_models(models);
assert(test_nj == nj)
assert(isequal(test, p))

% TODO Check ValueErrors occur
model_a = drss(2,1,1);


%% Check non-linear models

% Non-linear system model
model = struct();
model.n = 2;
model.nu = 0;
model.ny = 1;
model.f = @vdpStateFcn;
model.h = @vdpMeasurementFcn;
model.dfdx = @vdpStateJacobianFcn;
model.dhdx = @vdpMeasurementJacobianFcn;
model.Ts = 0.05;

[n, nu, ny, Ts] = check_model_NL(model);
assert(isequal([n, nu, ny], [2, 0, 1]))
assert(Ts == 0.05)

models = repmat({model}, 1, 3);
[nj, n, nu, ny, Ts] = check_models_NL(models);
assert(isequal([nj, n, nu, ny], [3, 2, 0, 1]))
assert(Ts == 0.05)