% plotting the smoothing function

clear all

addpath('../General_functions/')

t  = linspace(-2*pi, 2*pi, 1000);

C1 = 1;
C2 = pi;

f = smoothing_function(t, C1, C2);

figure()
plot(t, f)
title('Smoothing function', 'FontSize', 18)
xlabel('$t$', 'FontSize', 18)
ylabel('$f(t)$', 'FontSize', 18)
ylim([-0.1 1.1])