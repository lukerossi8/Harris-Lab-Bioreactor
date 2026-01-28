clearvars; close all; clc;
[X,Y] = meshgrid(0:0.2:2*pi); % X must be >= 0, Y must be <= 2pi
% Assume inviscid fluid for now
F_t = 1;
U = sin(X) .* cos(Y) * F_t;
V = -cos(X) .* sin(Y) * F_t;
[curlz,cav] = curl(X, Y, U, V);
c = pcolor(X, Y, curlz);
c.FaceColor = 'interp';
hold on
q = quiver(X, Y, U, V, 1.5, 'k');
xlabel('x position')
ylabel('y position')
title('Taylor Green Vortex Vector Field, Shaded With Vorticity')
grid on