% m1_m2_int test
% 

clear all
%
%
%
warning('PIM_hankel.m has been changed so that it containsm1 and m2 functions')
%
%
%

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

t = [a+h/2: h: b - h/2];
t_grid = [a:h:b];

err1 = zeros(size(s));
err2 = zeros(size(s));

f = @(x, y) 1i*besselh(0, k*abs(x - y))/4;

for j = 1:length(s)
    
    mat_approx(j) = integral(@(y) f(s(j), y), a, b);

    our_approx(j) = sum(m1(k, s(j), t, C1, C2).*w1_weights(k, s(j), ...
        t_grid(1:end - 1), t_grid(2:end))...
        + h*m2(k, s(j), t, C1, C2));

    PIM_fun = PIM_int_hankel_f(k, s(j), h, t, 1, t_grid, C1, C2);

    err1(j) = mat_approx(j) - our_approx(j);

    err2(j) = mat_approx(j) - PIM_fun;


end

figure()
plot(s, err1, 'DisplayName', 'Modular code', 'LineStyle', '--')
hold on
plot(s, err2, 'DisplayName', 'Not modular code', 'LineStyle', '-.')
legend show
title('Plot of error difference')

figure()
plot(s, mat_approx, 'DisplayName', 'Matlab approximation')
hold on 
plot(s, our_approx, 'DisplayName', 'Our approximation')
hold off
legend show
title('Comparison plot')

%% narrowing down the potential error region

clear s
s = linspace(a-0.0001, b+0.0001, 1000);

err_a_b = zeros(size(s));


for j = 1:length(s)
    
    mat_approx_a_b(j) = integral(@(y) f(s(j), y), a, b);

    our_approx_a_b(j) = sum(m1(k, s(j), t, C1, C2).*w1_weights(k, s(j), ...
        t_grid(1:end - 1), t_grid(2:end))...
        + h*m2(k, s(j), t, C1, C2));


    err_a_b(j) = mat_approx_a_b(j) - our_approx_a_b(j);


end

figure()
plot(s, err_a_b)
title('Differnce between Matlab and our approximation to  $\frac{i}{4} \int_{a}^{b} H_{0}^{(1)}(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('Difference', 'fontsize',15,'interpreter','latex')

figure()
plot(s, mat_approx_a_b, 'DisplayName', 'Matlab approximation')
hold on 
plot(s, our_approx_a_b, 'DisplayName', 'Our approximation')
hold off
legend show
title('Comparison of different approximations to $\frac{i}{4} \int_{a}^{b} H_{0}^{(1)}(k \vert s - t \vert) \mathrm{d} t$','fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('I(s)', 'fontsize',15,'interpreter','latex')


%% looking at a smaller subset still
a_test = s(435);
b_test = s(445);
s_sing_test = linspace(a_test, b_test, 1000);
err_sing_test = zeros(size(s_sing_test));

% t = [a_test+h/2:h:b_test- h/2];
% t_grid = [a_test:h:b_test];



for j = 1:length(s)
    
    mat_approx_sing_test(j) = integral(@(y) f(s_sing_test(j), y), a, b);

    our_approx_sing_test(j) = sum(m1(k, s_sing_test(j), t, C1, C2).*w1_weights(k, s_sing_test(j), ...
        t_grid(1:end - 1), t_grid(2:end))...
        + h*m2(k, s_sing_test(j), t, C1, C2));


    err_sing_test(j) = mat_approx_sing_test(j) - our_approx_sing_test(j);


end

figure()
plot(s_sing_test, err_sing_test)
title('Isolating a problem case showing differnce between Matlab and our approximation to  $\frac{i}{4} \int_{a}^{b} H_{0}^{(1)}(k \vert s - t \vert) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('Difference', 'fontsize',15,'interpreter','latex')

figure()
plot(s_sing_test, mat_approx_sing_test, 'DisplayName', 'Matlab approximation')
hold on 
plot(s_sing_test, our_approx_sing_test, 'DisplayName', 'Our approximation')
hold off
legend show
title('Isolating a problem case comparing different approximations to $\frac{i}{4} \int_{a}^{b} H_{0}^{(1)}(k \vert s - t \vert) \mathrm{d} t$','fontsize',15,'interpreter','latex')
xlabel('$s$', 'fontsize',15,'interpreter','latex')
ylabel('I(s)', 'fontsize',15,'interpreter','latex')
