% Test make_trans_funcs_ssd.m

rng(0)
addpath("../process-observers")


%% Test with SISO 1st order system

sys_rodin_step

[f,h,params] = make_trans_funcs_ssd(A,B,C,D);

assert(isequal(params.A, A));
assert(isequal(params.B, B));
assert(isequal(params.C, C));
assert(isequal(params.D, D));

[n, nu, ny] = check_dimensions(A, B, C, D);
xk = randn(n, 1);
uk = randn(nu, 1);
dt = 0.5;

xkp1 = f(xk, uk, dt, params);
assert(isequal(xkp1, A * xk + B * uk));

yk = h(xk, uk, dt, params);
assert(isequal(yk, C * xk + D * uk));


%% Test with random system

n = 4;
nu = 3;
ny = 2;

sys = drss(n,ny,nu);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

[f,h,params] = make_trans_funcs_ssd(A,B,C,D);

xk = randn(n, 1);
uk = randn(nu, 1);
dt = 2;

xkp1 = f(xk, uk, dt, params);
assert(isequal(xkp1, A * xk + B * uk));

yk = h(xk, uk, dt, params);
assert(isequal(yk, C * xk + D * uk));
