% Test function string2latex.m

%% Test with strings

assert(strcmp(string2latex('abc'), '$abc$'))
assert(strcmp(string2latex('$abc$'), '$abc$'))
assert(strcmp(string2latex('$100'), '$$100$'))
assert(strcmp(string2latex('100 $'), '$100 $$'))


%% Check with cell arrays of strings

labels1 = {'y_1(k)', 'y_2(k)', 'y_3(k)'};
labels2 = string2latex(labels1);
assert(isequal(labels2, {'$y_1(k)$', '$y_2(k)$', '$y_3(k)$'}))
assert(isequal(string2latex(labels2), {'$y_1(k)$', '$y_2(k)$', '$y_3(k)$'}))
