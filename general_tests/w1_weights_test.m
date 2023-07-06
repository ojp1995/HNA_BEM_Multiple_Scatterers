% in this function we will be testing that the integrator w1_weights.m
% works, by computing the integral \int_a^b ln(k|s - t|) dt, for a range of
% values of s. To compare this we will be using the matlab inbuilt
% integrator.

clear all

addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')

k = 0.85;
lambda = 2*pi/k;
a = 1.35;
b = 1 + lambda/10;

N = 2^20;
h = (b - a)/N;
C1 = 1;
C2 = 2*pi;

t = [a+h/2:h:b- h/2];
t_grid = [a:h:b];


s = linspace(a + h/2, b+h/2, 1000);


for j = 1:length(s)
    t_in_s_q1(j) = sum(s(j) == t_grid);
end

err = zeros(size(s));

f = @(x, y) log(k*abs(x - y));

for j = 1:length(s)

    matlab_val(j) = integral(@(y) f(s(j), y), a, b);

    our_approx(j) = sum(w1_weights(k, s(j), t_grid(1:end - 1), t_grid(2:end)));

    err(j) = matlab_val(j) - our_approx(j);

end

% err
% 
% [matlab_val.' our_approx.']

figure()
plot(s, err)
title('Difference between Matlab integrator and analytic evaluation of $\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('Difference', 'fontsize',15,'interpreter','latex')

figure()
plot(s, matlab_val, 'DisplayName', 'Matlabs approximation')
hold on
plot(s, our_approx, 'DisplayName', 'Our approximation')
legend show
title('Comparison between Matlab integrator and analytic evaluation of $\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel( '$s$',  'fontsize',15,'interpreter','latex')
ylabel('$\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$' ,  'fontsize',15,'interpreter','latex')

%%
% what if we look at a really short subset, specifically within the bounds
% a and b. Then  look even further in detail to being at a very short
% subset of the values.

clear s
s = linspace(a-0.0001, b+0.0001, 1000);

err_a_b = zeros(size(s));

for j = 1:length(s)
    t_in_s_q2(j) = sum(s(j) == t_grid);
end

for j = 1:length(s)

    matlab_val_a_b(j) = integral(@(y) f(s(j), y), a, b);

    our_approx_a_b(j) = sum(w1_weights(k, s(j), t_grid(1:end - 1), t_grid(2:end)));

    err_a_b(j) = matlab_val(j) - our_approx(j);

end

figure()
plot(s, matlab_val_a_b, 'DisplayName', 'Matlabs approximation')
hold on
plot(s, our_approx_a_b, 'DisplayName', 'Our approximation')
legend show
title('Comparison between Matlab integrator and analytic evaluation of $\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel( '$s$',  'fontsize',15,'interpreter','latex')
ylabel('$\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$' ,  'fontsize',15,'interpreter','latex')


figure()
plot(s, err_a_b)
title('Difference between Matlab integrator and analytic evaluation of $\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('Difference', 'fontsize',15,'interpreter','latex')

% looking at an even smaller subset, values (435:445)
figure()
plot(s(435:445), matlab_val_a_b(435:445), 'DisplayName', 'Matlabs approximation')
hold on
plot(s(435:445), our_approx_a_b(435:445), 'DisplayName', 'Our approximation')
legend show
title('Small subset of singular intergal - Comparison between analytic evaluation and Matlab integrator', 'fontsize',15,'interpreter','latex')
xlabel( '$s$',  'fontsize',15,'interpreter','latex')
ylabel('$\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$' ,  'fontsize',15,'interpreter','latex')


figure()
plot(s(435:445), err_a_b(435:445))
title('Singular part only - Plot of the difference between Matlab integrator and analytic evaluation', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('Difference', 'fontsize',15,'interpreter','latex')

%%
% How are we sure that we are right and this isn't correc, lets throw more
% points at it:
a_test = s(435);
b_test = s(445);
s_sing_test = linspace(a_test, b_test, 1000);
err_sing_test = zeros(size(s_sing_test));

t = [a_test+h/2:h:b_test- h/2];
t_grid = [a_test:h:b_test];

for j = 1:length(s_sing_test)
    t_in_s_q3(j) = sum(s_sing_test(j) == t_grid);
end

for j = 1:length(s_sing_test)
    
    matlab_sing_test(j) = integral(@(y) f(s_sing_test(j), y), a_test, b_test);

    our_approx_sing_test(j) = sum(w1_weights(k, s_sing_test(j), t_grid(1:end - 1), t_grid(2:end)));

    err_sing_test(j) = matlab_val(j) - our_approx(j);
end

figure()
plot(s_sing_test, matlab_sing_test, 'DisplayName', 'Matlabs approximation')
hold on
plot(s_sing_test, our_approx_sing_test, 'DisplayName', 'Our approximation')
legend show
title('Isolating a problem area and comparing between Matlab integrator and analytic evaluation of $\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel( '$s$',  'fontsize',15,'interpreter','latex')
ylabel('$\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$' ,  'fontsize',15,'interpreter','latex')

figure()
plot(s_sing_test, err_sing_test)
title('Isolating a problem area and showing the difference between Matlab integrator and analytic evaluation of $\int_{a}^{b}\ln(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('Difference', 'fontsize',15,'interpreter','latex')

disp([sum(t_in_s_q1), sum(t_in_s_q2), sum(t_in_s_q3)])