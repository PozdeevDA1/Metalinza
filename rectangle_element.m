function rectangle = rectangle_element(geometry, params, position, label)

x = position(1);
y = position(2);

rectangle_tag = sprintf('element_%i', label);
rectangle = geometry.feature.create(rectangle_tag,'Rectangle');
rectangle.set('size', [params.size_x, params.size_y]);
rectangle.set('pos', [x, y]);
rectangle.set('base', 'corner');
rectangle.set('createselection', 'on');
