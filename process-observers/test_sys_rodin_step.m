%% Test sys_rodin_step.m

sys_rodin_step

% Check models are identical
t_test = 0.5*(0:10)';
U_test = reshape(idinput(11*2), 11, 2);
[y_test, t_test] = lsim(Gpss, U_test, t_test);
[y_test2, t_test] = lsim(Gpd, U_test, t_test);
assert(all(abs(y_test - y_test2) < 1e-6))
