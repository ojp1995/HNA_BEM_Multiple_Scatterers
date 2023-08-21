% PIM graded rescalled and no-graded comparison

clear all
addpath('../General_functions/')

a = 0;
b = 2*pi;
L = b - a;
k = 10;

C1 = 1;
C2 = pi;



Lgrad = L*0.15;  % is this reasonable/may need to be tweaked!

s = L - Lgrad + 0.367 ;  % may change to being vector or a range of different points later

f_mat = @(t) (1i/4)*besselh(0, k*abs(s - t))./sqrt(t.*(L - t));

f = @(t) 1./sqrt(t.*(L - t));

int_mat= integral(@(t) f_mat(t), a, b);

alpha = linspace(1, 6, 6);
N_init = 160;
N_it_max = 10;
% h_init = (b - a)/N_init;
Q = N_init;

PIM_approx_standard = zeros(N_it_max, 1);
PIM_graded = zeros(N_it_max, length(alpha));
graded_PIM_abs_err = zeros(N_it_max, length(alpha));

%% computing 'true solution'
% h_true = h_init*2^-15;
h_true = 1e-6;
ts_grid_true = [a:h_true:b];
ts_mid_true = (ts_grid_true(2:end) + ts_grid_true(1:end-1))/2;
% ts_w_true = h_true;

% midpoint approx part
PIM_approx_standard_true = PIM_int_hankel_f(k, s, h_true, ts_mid_true, ...
    1, ts_grid_true, C1, C2);

for h_n = 1:N_it_max
    h = L/Q(h_n);
    ts_grid = [a:h:b];
    ts_mid = (ts_grid(2:end) + ts_grid(1:end-1))/2;
    ts_w = h;

    % midpoint approx part
    PIM_approx_standard(h_n) = PIM_int_hankel_f(k, s, h, ts_mid, 1,...
        ts_grid, C1, C2);

    PIM_abs_err(h_n) = abs(PIM_approx_standard(h_n) - ...
        PIM_approx_standard_true);

     for a_n = 1:length(alpha)

         % computign the grid and quadrature points/weights
         [t_grid_graded, t_mid_graded, w_graded, h(h_n, a_n)]...
            = get_graded_midpoint_half_interval(L, Lgrad, Q(h_n),...
            alpha(a_n));

         PIM_graded1 = graded_PIM_int_hankel_f(k, s, w_graded,...
             t_mid_graded, f(t_mid_graded), t_grid_graded, C1, C2 );

          PIM_graded2 = graded_PIM_int_hankel_f(k, L-s, w_graded,...
             t_mid_graded, f(t_mid_graded), t_grid_graded, C1, C2 );

          PIM_graded(h_n, a_n) = PIM_graded1 + PIM_graded2;

          graded_PIM_abs_err(h_n, a_n) = abs(PIM_approx_standard_true ...
              - PIM_graded(h_n, a_n));

     end 
     Q(h_n + 1) = Q(h_n)*2;

end

graded_PIM_abs_err, PIM_abs_err, abs(int_mat - PIM_graded), abs(int_mat - PIM_approx_standard)
%% EOC computations



for h_n = 1:N_it_max-1  % loop for h stepping

    EOC_standard(h_n) = log2(PIM_abs_err(h_n)/PIM_abs_err(h_n+1));

    for a_n = 1:length(alpha)

        EOC_graded(h_n, a_n) = log2(graded_PIM_abs_err(h_n, a_n)/...
            graded_PIM_abs_err(h_n+1, a_n));

    end
end

EOC_standard, EOC_graded