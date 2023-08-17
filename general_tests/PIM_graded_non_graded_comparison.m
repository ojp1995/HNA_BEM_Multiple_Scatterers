% PIM with graded midpoint test
% We will be initially running three comparisons. Firstly the current
% midpoint test, then the graded routine and then the matlab routine to the
% function f(s) = (i/4)int_0^L H_0^(1)(k|s - t|) dt. Then we will introduce
% the extra function (t(L - t))^{-1/2} into the kernel to mimic the phi
% function.
%
% We also want convergence results

clear all

addpath('../General_functions/')

a = 0;
b = 2*pi;
L = b - a;
k = 10;

C1 = 1;
C2 = pi;



Lgrad = L*0.15;  % is this reasonable/may need to be tweaked!

s = L - Lgrad - 0.367 ;  % may change to being vector or a range of different points later

f = @(t) (1i/4)*besselh(0, k*abs(s - t))./sqrt(t.*(L - t));

int_mat= integral(@(t) f(t), a, b);

% int_WA = -0.12382106204002354372091606659770831438124999253547176499371697401508783434823130097297493804577529998680421785666968663001...
%     +1i*0.67376315733699508561573656398060748211027002613870528634210843805787874105337575848695338655363343713758030633756225240020;

% int_WA = -0.12382106204002354372091606659770831438124999253547176499371697...
%     + 1i*0.67376315733699508561573656398060748211027002613870528634210843;

int_WA_k10_singularites = 1i*(0.0626194 - 0.0636219*1i)/4;
alpha = linspace(1, 6, 6);
N_init = 160;
N_it_max = 10;
h_init = (b - a)/N_init;

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

%% 
h =zeros(N_it_max, 1);
h(1) = h_init;
for h_n = 1:N_it_max  % loop for h stepping
    % all info for standard midpoint approximation
    ts_grid = [a:h(h_n):b];
    ts_mid = (ts_grid(2:end) + ts_grid(1:end-1))/2;
    ts_w = h(h_n);

    % midpoint approx part
    PIM_approx_standard(h_n) = PIM_int_hankel_f(k, s, h(h_n), ts_mid, 1,...
        ts_grid, C1, C2);

    PIM_abs_err(h_n) = abs(PIM_approx_standard(h_n) - ...
        PIM_approx_standard_true);

    for a_n = 1:length(alpha)

        [t_grid_graded1, t_mid_graded1, w_graded1, ~, t_grid_graded2, ...
            t_mid_graded2, w_graded2] = ...
            get_graded_midpoint_quad_points_split(L, Lgrad, ...
            h(h_n), alpha(a_n));

        PIM_graded1 = graded_PIM_int_hankel_f(k, s, w_graded1,...
            t_mid_graded1, 1, t_grid_graded1, C1, C2);

        % I believe something is wrong in here! Test to make sure that
        % there is agreement with the above integral! That will let us know
        % that there is a substitution error here somewhere.
%         if s < L - Lgrad
            PIM_graded2 = graded_PIM_int_hankel_f(k, L - s, w_graded2,...
                t_mid_graded2, 1, t_grid_graded2, C1, C2);
%         else
%             PIM_graded2 = graded_rescalled_PIM_int_hankel_f(k, s, ...
%                 w_graded2, t_mid_graded2, 1, t_grid_graded2, C1, C2, L);
% 
% 
%         end

        PIM_graded(h_n, a_n) = PIM_graded1 + PIM_graded2;

        graded_PIM_abs_err(h_n, a_n) = abs(PIM_graded(h_n, a_n)...
            - PIM_approx_standard_true);

    end 
    h(h_n+1) = h(h_n)/2;

end
graded_PIM_abs_err, PIM_abs_err %, abs(int_mat - PIM_graded), abs(int_mat - PIM_approx_standard)
%% EOC computations



for h_n = 1:N_it_max-1  % loop for h stepping

    EOC_standard(h_n) = log2(PIM_abs_err(h_n)/PIM_abs_err(h_n+1));

    for a_n = 1:length(alpha)

        EOC_graded(h_n, a_n) = log2(graded_PIM_abs_err(h_n, a_n)/...
            graded_PIM_abs_err(h_n+1, a_n));

    end
end

EOC_standard, EOC_graded