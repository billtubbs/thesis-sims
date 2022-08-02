% Test matrix_element_labels.m

clear

%% Test 1

labels = matrix_element_labels('a', 1:2, 1:3, '');
labels1 = {'a_1_1', 'a_1_2', 'a_1_3';
      'a_2_1', 'a_2_2', 'a_2_3'};
assert(isequal(labels, labels1))


%% Test 2

labels = matrix_element_labels('x', {'a','b','c'}, {1,2}, '(t)');
labels2 = {'x_a_1(t)', 'x_a_2(t)';
           'x_b_1(t)', 'x_b_2(t)';
           'x_c_1(t)', 'x_c_2(t)'};
assert(isequal(labels, labels2))


%% Test 3

labels = matrix_element_labels('\hat{x}', {'a','b','c'}, {1,2}, '(t)', true);
labels3 = {'$\hat{x}_a_1(t)$', '$\hat{x}_a_2(t)$';
           '$\hat{x}_b_1(t)$', '$\hat{x}_b_2(t)$';
           '$\hat{x}_c_1(t)$', '$\hat{x}_c_2(t)$'};
assert(isequal(labels, labels3))
