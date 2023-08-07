% graded midpoint LN SINGULARITY test
% We are going to compute the integral int_0^2 ln(t)(t(2 - t))^-1/2 dt. The
% answer is taken to be -pilog(2) https://www.wolframalpha.com/input?i2d=true&i=Integrate%5BDivide%5Bln%5C%2840%29t%5C%2841%29%2CSqrt%5Bt%5C%2840%292-t%5C%2841%29%5D%5D%2C%7Bt%2C0%2C2%7D%5D


clear all
addpath('../General_functions/')

a = 0;
b = 2;
L = b-a;
Lgrad = L*0.15;  % is this reasonable/may need to be tweaked!


int_true = -pi*log(2);

alpha = linspace(1, 10, 10);
N_init = 160;
N_it_max = 13;
h = (b - a)/N_init;

f = @(t) log(abs(1 - t))./sqrt(t.*(2 - t));

for h_n = 1:N_it_max  % loop for h stepping

    for a_n = 1:length(alpha)
        [t_grid, t_mid, w, Q(h_n, a_n)] = ...
            get_graded_midpoint_quad_points(L, Lgrad, h(h_n), alpha(a_n));

        [t_grid_graded1, t_mid_graded1, w_graded1, ~, t_grid_graded2, ...
            t_mid_graded2, w_graded2] = ...
            get_graded_midpoint_quad_points_split(L, Lgrad, ...
            h(h_n), alpha(a_n));

        mid_point_approx(h_n, a_n) = sum(w.*f(t_mid));

        midpoint_f_split_approx(h_n, a_n) = ...
            sum(w_graded1.*f(t_mid_graded1)) ...
            + sum(w_graded2.*f(t_mid_graded2));

        err_abs(h_n, a_n) =  abs(mid_point_approx(h_n, a_n) - int_true);

        err_abs_graded(h_n, a_n) =  ...
            abs(midpoint_f_split_approx(h_n, a_n) - int_true);

        err_abs_norm(h_n, a_n) =  abs(mid_point_approx(h_n, a_n) ...
            - int_true)/abs(int_true);

        clear t_grid t_mid

    end
    h(h_n+1) = h(h_n)/2;
end

for h_n = 1:N_it_max-1  % loop for h stepping

    for a_n = 1:length(alpha)
        
        EOC(h_n, a_n) = log2(err_abs(h_n, a_n)...
            /err_abs(h_n + 1, a_n));

        EOC_graded(h_n, a_n) = log2(err_abs_graded(h_n, a_n)...
            /err_abs_graded(h_n + 1, a_n));
    end 
end

err_abs, err_abs_norm, Q, EOC, h


err_abs_graded, EOC_graded
