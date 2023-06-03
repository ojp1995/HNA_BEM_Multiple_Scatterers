% script to test that m1 and m2 are working as they should be
%
% We will fix s at a range of values and then comoute a range of values for
% t, between 0 and L and compare to make sure that:
% i*H_0^(1)(k|s - t|)/4 = sigma1(s, t)m1(s, t) + m2(s, t)
clear all


addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
k = 0.85;
lambda = 2*pi/k;
a = 1.35;
b = 1 + lambda/10;

s = linspace(0, 10, 1000);

N = 2^10;
h = (b - a)/N;
C1 = 1;
C2 = 2*pi;

t = a + 0.3;



err = zeros(size(s));

% f = @(x, y) 1i*besselh()
for j = 1:length(s)
    dist = abs(s(j) - t);

    hankel(j) = 1i*besselh(0, k*dist)/4;
    our_approx(j) = m1(k, s(j), t, C1, C2).*log(k*abs(s(j) - t))...
        + m2(k, s(j), t, C1, C2);

    err(j) = hankel(j) - our_approx(j);

end

figure()
plot(s, err)
title('Difference between Matlabs Bessel \\ function and deconstruction','fontsize',15,'interpreter','latex')
xlim([-0.1, 10.1])
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('$\frac{i}{4}H_{0}^{(1)}(k \vert s - t \vert) - (m_{1}(s, t)\sigma_{1}(s, t) + m_{2}(s, t) \sigma_{2}(s, t))$', 'fontsize',15,'interpreter','latex')

figure()
plot(s, hankel, 'DisplayName', 'Matlab approximation', 'LineStyle', '--')
hold on
plot(s, our_approx, 'DisplayName', 'Deconstruction', 'LineStyle', '-.')
hold off
legend show
title('Comparison between our approx of bessel function and Matlabs', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('$f(s, t)$' , 'fontsize',15,'interpreter','latex')