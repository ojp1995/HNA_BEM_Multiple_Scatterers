% In this function we will be showing the poor performance of the code when
% increasing the wavenumber

clear all

% only adding one path to centralise solvers
addpath('../General_functions/')

G1_data.G = [-2*pi, 2*pi, 0, 0];

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 

Lgrad_coeff = 0.15;
alpha = 4;

C_wl= 1/20;

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;


k = [10, 50, 100];

N_bf = [10, 20, 40, 80, 160];

for kk = 1:length(k)
    [G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N, G1_data.L] = ...
        discretistion_vars_graded(G1_data.G, C_wl, k(kk), Lgrad_coeff, alpha);

    [G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, G2_data.t_grid, G2_data.t_mid, G2_data.w, G2_data.N, G2_data.L] = ...
        discretistion_vars_graded(G2_data.G, C_wl, k(kk), Lgrad_coeff, alpha);


    for n = 1:length(N_bf)
        tic
        [x1_1, y1_1, x1_2, y1_2, t1_bf_grid, t1_mid, ~, ~] = ...
        given_N_discretistion_vars_graded(G1_data.G, N_bf(n), Lgrad_coeff, alpha);
        [x2_1, y2_1, x2_2, y2_2, t2_bf_grid, t2_mid, ~, ~] = ...
        given_N_discretistion_vars_graded(G2_data.G, N_bf(n), Lgrad_coeff, alpha);
    
        % collocation points
        G1_data.s = [ t1_mid ; flip(G1_data.L - t1_mid) ];
        G1_data.x_col = [ x1_1 ; flip(x1_2) ];
        G1_data.y_col = [ y1_1 ; flip(y1_2) ];
        
        % col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
        G2_data.s = [ t2_mid ; flip(G2_data.L - t2_mid) ];
        G2_data.x_col = [ x2_1 ; flip(x2_2) ];
        G2_data.y_col = [ y2_1 ; flip(y2_2) ];

        [S11, S12, S21, S22, u_inc1, u_inc2] = ...
        compute_matrices_for_iterative_solve(G1_data, G2_data, k(kk), ...
        t1_bf_grid, t2_bf_grid, theta, C1, C2 );
    
        % all in one solve
        A = [S11  S12 ; S21 S22];
        
        u_inc = [u_inc1 ; u_inc2];
        
        coeffs = A\u_inc;
    
        t1_bf_grid_store{kk, n, :} = t1_bf_grid;
        t2_bf_grid_store{kk, n, :} = t2_bf_grid;
        
        aj_1{kk, n, :} = coeffs(1:2*length(t1_bf_grid) - 2);
        aj_2{kk, n, :} = coeffs(2*length(t1_bf_grid) - 2+ 1: end);

        clear A S11 S12 S21 S22 uinc u_inc1 u_inc2

        toc
    



    end 
    
end


%% plotting



for kk = 1:length(k)
    C_wl_err = 1/40;

    [x1_1_q_err, y1_1_q_err, x1_2_q_err, y1_2_q_err, t1_grid_err, ...
        t1_mid_err, w1_err, N1_err, L1] = ...
        discretistion_vars_graded(G1_data.G, C_wl_err, k(kk), Lgrad_coeff, alpha);
    
    [x2_1_q_err, y2_1_q_err, x2_2_q_err, y2_2_q_err, t2_grid_err, ...
        t2_mid_err, w2_err, N2_err, L2] = ...
        discretistion_vars_graded(G2_data.G, C_wl_err, k(kk), Lgrad_coeff, alpha);
    
    t1_mid_err = t1_mid_err(100:end);
    t2_mid_err = t2_mid_err(100:end);
    
    x1_plotting = [t1_mid_err; flip(L1 - t1_mid_err) ];
    x2_plotting = [t2_mid_err; flip(L2 - t2_mid_err) ];
    
    w1_err = [w1_err(100:end); flip(w1_err(100:end))];
    w2_err = [w2_err(100:end); flip(w2_err(100:end))];
    
%     phi_1 = zeros(length(N_bf), length(t1_mid_err)*2);
%     phi_2 = zeros(length(N_bf), length(t2_mid_err)*2);

    for j = 1:length(N_bf)

        phi_1{kk, j, :} = graded_coeff_2_solution(aj_1{kk, j}, t1_bf_grid_store{kk, j},...
            t1_mid_err, G1_data.L);
        phi_2{kk, j, :} =graded_coeff_2_solution(aj_2{kk, j}, t2_bf_grid_store{kk, j},...
            t2_mid_err, G2_data.L);
    end

     figure()
    for j = 1:length(N_bf)
        plot(x1_plotting/G1_data.L, phi_1{kk, j}, 'DisplayName', ...
            strcat('N = ', num2str(length(aj_1{kk, j}))))
        hold on

    end
    xlabel('$x/L_1$')
    ylabel('$\phi_1(x)$')
    title(['Polynoial approximation solve, $\phi_{1}(x)$, k = ', num2str(k(kk))])
    legend show

    figure()
    for j = 1:length(N_bf)
        plot(x2_plotting/G2_data.L, phi_2{kk, j}, 'DisplayName', ...
            strcat('N = ', num2str(length(aj_2{kk, j}))))
        hold on

    end
    xlabel('$x/L_2$')
    ylabel('$\phi_2(x)$')
    title(['Polynoial approximation solve, $\phi_{2}(x)$, k = ', num2str(k(kk))])
    legend show

    % error and convergence
    phi_1_only_norm = sum(w1_err.*abs(phi_1{kk, end}));
    phi_2_only_norm = sum(w2_err.*abs(phi_2{kk, end}));
    
    
    for j = 1:length(N_bf) - 1
        
        err_phi_1(kk, j) = sum(w1_err.*abs(phi_1{kk, end} - phi_1{kk, j}))./phi_1_only_norm;
            
    
        err_phi_2(kk, j) = sum(w2_err.*abs(phi_2{kk, end} - phi_2{kk, j}))./phi_2_only_norm;
    
    end


end
% %
% for kk = 1:length(k)
%     figure()
%     for j = 1:length(N_bf)
%         plot(x1_plotting/G1_data.L, phi_1{kk, j}, 'DisplayName', ...
%             strcat('N = ', num2str(length(aj_1{kk, j}))))
%         hold on
% 
%     end
%     xlabel('$x/L_1$')
%     ylabel('$\phi_1(x)$')
%     title(['Polynoial approximation solve, $\phi_{1}(x)$, k = ', num2str(kk)])
%     legend show
% 
% end
% 
% for kk = 1:length(k)
%     figure()
%     for j = 1:length(N_bf)
%         plot(x2_plotting/G2_data.L, phi_2{kk, j}, 'DisplayName', ...
%             strcat('N = ', num2str(length(aj_2{kk, j}))))
%         hold on
% 
%     end
%     xlabel('$x/L_2$')
%     ylabel('$\phi_2(x)$')
%     title(['Polynoial approximation solve, $\phi_{2}(x)$, k = ', num2str(kk)])
%     legend show
% 
% end

%% Error/ convergence plot

% for kk = 1:length(k)
%     phi_1_only_norm = sum(w1_err.'.*abs(phi_1{kk, end}));
%     phi_2_only_norm = sum(w2_err.'.*abs(phi_2{kk, end}));
%     
%     
%     for j = 1:length(N_bf) - 1
%         
%         err_phi_1(kk, j) = sum(w1_err.'.*abs(phi_1{kk, end} - phi_1{kk, j}))./phi_1_only_norm;
%             
%     
%         err_phi_2(kk, j) = sum(w2_err.'.*abs(phi2{kk, end} - phi_1{kk, j}))./phi_2_only_norm;
%     
%     end
% 
% end

%% EOC

for kk = 1:length(k)

    for j = 1:length(N_bf) - 2

        EOC_phi_1(kk, j) = log2(err_phi_1(kk, j)/err_phi_1(kk, j+1));
    
        EOC_phi_2(kk, j) = log2(err_phi_2(kk, j)/err_phi_2(kk, j+1));
    end

end

err_phi_1, err_phi_2, EOC_phi_1, EOC_phi_2


