%graded midpoint integral test (new quadrature)
% We are going to compute the integral int_0^2 (t(2 - t))^-1/2 dt. The
% answer we are taking to be true is pi (https://www.wolframalpha.com/input?i2d=true&i=Integrate%5BPower%5B%5C%2840%29t%5C%2840%292+-+t%5C%2841%29%5C%2841%29%2C-Divide%5B1%2C2%5D%5D%2C%7Bt%2C0%2C2%7D%5D)

clear all
addpath('../General_functions/')

a = 0;
b = 2;
L = b-a;
Lgrad = L*0.15;  % is this reasonable/may need to be tweaked!

int_true = pi;

alpha = linspace(1, 10, 10);
N_init = 160;
N_it_max = 15;
% h = (b - a)/N_init;
Q = N_init;

f = @(t) 1./sqrt(t.*(L - t));

matlab_int = integral(f, 0, L);
for h_n = 1:N_it_max  % loop for h stepping

    for a_n = 1:length(alpha)
        [t_grid_gradedLHS, t_mid_gradedLHS, w_graded_LHS, hLHS(h_n, a_n)]...
            = get_graded_midpoint_half_interval(L, Lgrad, Q(h_n),...
            alpha(a_n));

        [t_grid_gradedRHS, t_mid_gradedRHS, w_graded_RHS, hRHS(h_n, a_n)] =...
            get_graded_midpoint_secondhalf_interval(L, Lgrad, Q(h_n), ...
            alpha(a_n));
        midpoint_graded(h_n, a_n) = sum( w_graded_LHS.*f(t_mid_gradedLHS) ...
            + w_graded_RHS.*f(t_mid_gradedRHS));

        midpoint_graded_rescalled(h_n, a_n) = sum( ...
            w_graded_LHS.*f(t_mid_gradedLHS) ...
            + w_graded_LHS.*f(t_mid_gradedLHS));
%             + flipud(w_graded_LHS).*f(L/2 + flipud(t_mid_gradedLHS)));

        err_graded(h_n, a_n) = abs(matlab_int - midpoint_graded(h_n, a_n));

        err_graded_rescalled(h_n, a_n) = abs(matlab_int - ...
            midpoint_graded_rescalled(h_n, a_n));

    end
    Q(h_n+1) = Q(h_n)*2;
end

% EOC computation
for h_n = 1:N_it_max-1
    for a_n = 1:length(alpha)
        EOC(h_n, a_n) = log2(err_graded(h_n, a_n)...
            /err_graded(h_n + 1, a_n));

        EOC_rescalled(h_n, a_n) = log2(err_graded_rescalled(h_n, a_n)...
            /err_graded_rescalled(h_n + 1, a_n));

    end 
end

err_graded, EOC, err_graded_rescalled, EOC_rescalled