% in this function we will be testing that the integrator w1_weights.m
% works, by computing the integral \int_a^b ln(k|s - t|) dt, for a range of
% values of s. To compare this we will be using the matlab inbuilt
% integrator.

clear all

% addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')

k = 0.85;
lambda = 2*pi/k;
a = 1.35;
b = 1 + lambda/10;

s = linspace(a - 0.2, b+0.2, 1000);

N = 2^20;
h = (b - a)/N;
C1 = 1;
C2 = 2*pi;

t = [a+h/2:h:b- h/2];
t_grid = [a:h:b];


err = zeros(size(s));

f = @(x, y) log(k*abs(x - y));

for j = 1:length(s)

    matlab_val(j) = integral(@(y) f(s(j), y), a, b);

    our_approx(j) = sum(w1_weights(k, s(j), t_grid(1:end - 1), t_grid(2:end)));

    err(j) = sum(matlab_val(j) - our_approx(j));

end

% err
% 
% [matlab_val.' our_approx.']

figure()
plot(s, err)
title('Plot of the difference between Matlabs integrator and our approximation')

figure()
plot(s, log(matlab_val), 'DisplayName', 'Matlabs approximation')
hold on
plot(s, log(our_approx), 'DisplayName', 'Our approximation')
legend show
title('Comparison between our approximation and matlabs')

%%
% what if we look at a really short subset, specifically within the bounds
% a and b. Then  look even further in detail to being at a very short
% subset of the values.

clear s
s = linspace(a-0.0001, b+0.0001, 1000);

err_a_b = zeros(size(s));


for j = 1:length(s)

    matlab_val_a_b(j) = integral(@(y) f(s(j), y), a, b);

    our_approx_a_b(j) = sum(w1_weights(k, s(j), t_grid(1:end - 1), t_grid(2:end)));

    err_a_b(j) = sum(matlab_val(j) - our_approx(j));

end

figure()
plot(s, (matlab_val_a_b), 'DisplayName', 'Matlabs approximation')
hold on
plot(s, (our_approx_a_b), 'DisplayName', 'Our approximation')
legend show
title('Singular part only - Comparison between our approximation and matlabs')


figure()
plot(s, err_a_b)
title('Singgular part only - Plot of the difference between Matlabs integrator and our approximation')


% looking at an even smaller subset, values (435:445)
figure()
plot(s(435:445), (matlab_val_a_b(435:445)), 'DisplayName', 'Matlabs approximation')
hold on
plot(s(435:445), (our_approx_a_b(435:445)), 'DisplayName', 'Our approximation')
legend show
title('Small subset of singular intergal - Comparison between our approximation and matlabs')


figure()
plot(s, err_a_b)
title('Singgular part only - Plot of the difference between Matlabs integrator and our approximation')
