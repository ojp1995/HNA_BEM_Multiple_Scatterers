% Midpoint, different screens, comparison between graded quadrtature and
% non-graded quadrature

clear all

addpath('../General_functions/')

G1 = [-2*pi, 2*pi, 0, 0];

G2 = [2*pi, 0, 5*pi, 3*pi];

Lgrad_coeff = 0.15;
% h = 0.2;
k = 10;
Q = 100;
alpha = 2;
C_wl = 1/5; 

fnq = @(t, L) 1./sqrt(t.*(L - t));

% Quadratuer nodes for half the interval
[x1_1, y1_1, x1_2, x1_2, t1, t1_mid, w1, N1, L1] = ...
    discretistion_vars_graded( G1, C_wl, k, Lgrad_coeff, alpha);

[x2_1, y2_1, x2_2, y2_2, t2, t2_mid, w2, N2, L2] = ...
    discretistion_vars_graded( G2, C_wl, k, Lgrad_coeff, alpha);

% keyboard

% computing the graded_midpoint approximation
col_point = 7;

I_graded_approx = midpoint_hankel_f_diff_screen(k, x1_1,...
    y1_1, x2_1, y2_1, w2, fnq(t2_mid, L2)) + ...
    midpoint_hankel_f_diff_screen(k, x1_1, y1_1, x2_2,...
    y2_2, w2, fnq(t2_mid, L2));

[y1t, y2t, ~, t_mid, h, ~, ~, ~] = discretisation_variables(G2, C_wl/100, k);
I_midpoint = midpoint_hankel_f_diff_screen(k, x1_1, ...
    y1_1, y1t, y2t, h, fnq(t_mid, L2));

[I_graded_approx, I_midpoint]

sum(I_graded_approx - I_midpoint)/length(I_midpoint)

% Now we are sure it works we should test to see how changing alpha and
% C_wl effects the behaviour and also we should see at which point it
% breaks due to handling the small numbers not near zero.

alpha = 1:6;
C_wl = 1/5;
N_it = 10;  % number of iterations we will be looping over

% true
[y1t, y2t, ~, t_mid, h, ~, N_true, ~] = discretisation_variables(G2, C_wl/1000, k);
I_true = midpoint_hankel_f_diff_screen(k, x1_1, ...
    y1_1, y1t, y2t, h, fnq(t_mid, L2));
%%
% true graded
[x2_1, y2_1, x2_2, y2_2, t2, t2_mid, w2, N_ref, L2] = ...
    discretistion_vars_graded( G2, C_wl/(2^(N_it+1)), k, Lgrad_coeff, 7);
I_graded_true = midpoint_hankel_f_diff_screen(k, x1_1,...
    y1_1, x2_1, y2_1, w2, fnq(t2_mid, L2)) + ...
    midpoint_hankel_f_diff_screen(k, x1_1, y1_1, x2_2,...
    y2_2, w2, fnq(t2_mid, L2));
%%
for n = 1:N_it
    for a = 1:length(alpha)
        [x2_1, y2_1, x2_2, y2_2, t2, t2_mid, w2, N2(n, a), L2] = ...
        discretistion_vars_graded( G2, C_wl, k, Lgrad_coeff, alpha(a));

        I_graded_approx = midpoint_hankel_f_diff_screen(k, x1_1,...
        y1_1, x2_1, y2_1, w2, fnq(t2_mid, L2)) + ...
        midpoint_hankel_f_diff_screen(k, x1_1, y1_1, x2_2,...
        y2_2, w2, fnq(t2_mid, L2));

        I_graded_sum(n, a) = sum(I_graded_approx);

        err_graded(n , a) = sum(abs(I_true - I_graded_approx))/length(I_true);

        err_graded_graded_true(n, a) = sum(abs(I_graded_true - I_graded_approx))/length(I_graded_true);

    end
    C_wl = C_wl/2;
end


% estimate order of convergence

for n = 1:N_it-1
    for a = 1:length(alpha)

        EOC_graded(n , a) = log2(err_graded(n , a)/err_graded(n+1 , a) );
        EOC_graded_graded_true(n , a) = log2(err_graded_graded_true(n , a)/err_graded_graded_true(n+1 , a) );

    end
end

err_graded
EOC_graded
err_graded_graded_true
EOC_graded_graded_true


