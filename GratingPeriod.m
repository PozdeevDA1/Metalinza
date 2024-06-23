function [p, th] = GratingPeriod(params)

freq = params.wave.freq0;
c_w = params.wave.c_water;
lamb = c_w/freq;
focus = params.metalens.focus;
N = params.metalens.N;
x0 = params.metalens.x0;

x = zeros(1, N);
p = zeros(1, N);
th = zeros(1, N);


x(1) = x0;

syms p1
eqn = p1^4 + 2*x(1)*p1^3 + p1^2*(x(1)^2 - lamb^2) - 2*lamb^2*x(1)*p1 - lamb^2*(x(1)^2 + focus^2);
b =  vpasolve(eqn == 0, p1);
p(1) = b(imag(b) == 0 & real(b) > 0);
th(1) = atan((x(1) + p(1))/focus);

for n = 2:1:N
  x(n) = x(n-1) + p(n-1);
  syms pn
  eqn = pn^4 + 2*x(n)*pn^3 + pn^2*(x(n)^2 - lamb^2) - 2*lamb^2*x(n)*pn - lamb^2*(x(n)^2 + focus^2);
  b =  vpasolve(eqn == 0, pn);
  p(n) = b(imag(b) == 0 & real(b) > 0);
  th(n) = atan((x(1) + sum(p))/focus);
end
% 
% for idx_p = 1:length(p)
%     fprintf('n = %i, p = %.5f mm, th = %.2f \n', idx_p, p(idx_p)*1e3, rad2deg(th(idx_p)));
% end
% 
% 
% plot(1:N, rad2deg(th))
% xlabel('n', 'interpreter', 'latex');
% ylabel('$\theta_n$', 'interpreter', 'latex')
