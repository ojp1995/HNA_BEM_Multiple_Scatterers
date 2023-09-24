% polynomial all in one convergence test

clear all

% only adding one path to centralise solvers
addpath('../General_functions/')

% switch and things get interesting
G1_data.G = [-2*pi, 2*pi, 0, 0];
% G1_data.L = sqrt( (G1(3) - G1(1))^2 +(G1(4) - G1(2))^2 );

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 
% G2_data.L = sqrt( (G2(3) - G2(1))^2 +(G2(4) - G2(2))^2 );

Lgrad_coeff = 0.15;
alpha = 4;

C_wl= 1/10;

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

run_max = 5;

% quadrature nodes and other information needed
[G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N, G1_data.L] = ...
    discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);

[G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, G2_data.t_grid, G2_data.t_mid, G2_data.w, G2_data.N, G2_data.L] = ...
    discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);

G1_data.s = [ G1_data.t_mid(1:end) ; flip(G1_data.L - G1_data.t_mid(1:end)) ];
G1_data.x_col = [ G1_data.x_1_q(1:end) ; flip(G1_data.x_2_q(1:end)) ];
G1_data.y_col = [ G1_data.y_1_q(1:end) ; flip(G1_data.y_2_q(1:end)) ];

% col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
G2_data.s = [ G2_data.t_mid(1:end) ; flip(G2_data.L - G2_data.t_mid(1:end)) ];
G2_data.x_col = [ G2_data.x_1_q(1:end) ; flip(G2_data.x_2_q(1:end)) ];
G2_data.y_col = [ G2_data.y_1_q(1:end) ; flip(G2_data.y_2_q(1:end)) ];


for j = [1:5]
    % The support for the basis functions
    C_wl_bf1 = 1/(j+1);
    C_wl_bf2 = 1/(j+1);
    [~, ~, ~, ~, t1_bf_grid, ~, ~, ~, ~] = discretistion_vars_graded(...
        G1_data.G, C_wl_bf1, k, Lgrad_coeff, alpha);
    [~, ~, ~, ~, t2_bf_grid, ~, ~, ~, ~] = discretistion_vars_graded(...
        G2_data.G, C_wl_bf2, k, Lgrad_coeff, alpha);
    
    [S11, S12, S21, S22, u_inc1, u_inc2] = ...
        compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
        t1_bf_grid, t2_bf_grid, theta, C1, C2 );
    
    % all in one solve
    A = [S11  S12 ; S21 S22];
    
    u_inc = [u_inc1 ; u_inc2];
    
    coeffs = A\u_inc;

    t1_bf_grid_store{j, :} = t1_bf_grid;
    t2_bf_grid_store{j, :} = t2_bf_grid;
    
    aj_1{j, :} = coeffs(1:2*length(t1_bf_grid) - 2);
    aj_2{j, :} = coeffs(2*length(t1_bf_grid) - 2+ 1: end);

end


%%
x1_plot = linspace(0.01, G1_data.L/2 - 0.01, 40);

x1_plotting = [x1_plot(:) ; (G1_data.L/2 + x1_plot(:) )];

x2_plot = linspace(0.01, G2_data.L/2 - 0.01, 40);

x2_plotting = [x2_plot(:) ; (G2_data.L/2 + x2_plot(:) )];


phi_1 = zeros(run_max, length(x1_plotting));
phi_2 = zeros(run_max, length(x2_plotting));

for j = 1:5

    phi_1(j, :) = graded_coeff_2_solution(aj_1{j, :}, t1_bf_grid_store{j, :},...
        x1_plot, G1_data.L);
    phi_2(j, :) =graded_coeff_2_solution(aj_2{j, :}, t2_bf_grid_store{j, :},...
        x2_plot, G2_data.L);

end

% plotting

figure()
for j = 1:5
    plot(x1_plotting/G1_data.L, phi_1(j,:), 'DisplayName', ...
        strcat('C_wl = 1/', num2str(j+1)))
    hold on

end
xlabel('$x/L_1$')
ylabel('$\phi_1(x)$')
title('Polynoial approximation solve, $\phi_{1}(x)$')
legend show

figure()
for j = 1:5
    plot(x2_plotting/G2_data.L, phi_2(j,:), 'DisplayName', ...
        strcat('C_wl = 1/', num2str(j+1)))
    hold on

end
xlabel('$x/L_2$')
ylabel('$\phi_{2}(x)$')
title('Polynoial approximation solve, $\phi_{1}(x)$')
legend show

%% Error/ convergence plot

for j = 1:4

    err_phi_1(j) = sum(abs(phi_1(end, :) - phi_1(j, :))./abs(phi_1(end, :)))...
        /length(phi_1(end, :));

    err_phi_2(j) = sum(abs(phi_2(end, :) - phi_2(j, :))./abs(phi_2(end, :)))...
        /length(phi_2(end, :));

end

%% EOC

for j = 1:3

    EOC_phi_1(j) = log2(err_phi_1(j)/err_phi_1(j+1));

    EOC_phi_2(j) = log2(err_phi_2(j)/err_phi_2(j+1));

end

err_phi_1, err_phi_2, EOC_phi_1, EOC_phi_2