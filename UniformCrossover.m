function [y1, y2] = UniformCrossover(x1, x2, gamma0)

% Binary uniform crossover
% alpha = randi([0, 1], size(x1));

alpha = unifrnd(-gamma0, 1 + gamma0, size(x1));
y1 = alpha.*x1 + (1 - alpha).*x2;
y2 = alpha.*x2 + (1 - alpha).*x1;