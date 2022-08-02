% Test function vector_element_labels.m

label = vector_element_labels('y', '(k)', 1);
assert(isequal(label, {'y(k)'}))

label = vector_element_labels('y', '(k)', 1, true);
assert(isequal(label, {'$y(k)$'}))

labels = vector_element_labels('y', '(k)', 3);
assert(isequal(labels, {'y_1(k)', 'y_2(k)', 'y_3(k)'}))

labels = vector_element_labels('y', '(k)', 3, true);
assert(isequal(labels, {'$y_1(k)$', '$y_2(k)$', '$y_3(k)$'}))
