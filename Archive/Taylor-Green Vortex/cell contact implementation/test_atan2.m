a = [1, 0, 3, 4]; % [x1, x2, y1, y2]

dist_x = a(2) - a(1);
dist_y = a(4) - a(3);

theta = atan2(dist_y, dist_x);
theta
