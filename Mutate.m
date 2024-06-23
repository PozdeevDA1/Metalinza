function y = Mutate(x, mutation_rate, sigma0)

flag = rand(size(x)) < mutation_rate;

% y = x;
% y(flag) = 1 - x(flag);

%
r = randn(size(x));
step = sigma0 * r(flag);
y = x;
y(flag) = x(flag) + step;
%



end