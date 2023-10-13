% k -> \infty test idea

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

C_wl= 1/20;



theta = 0;

d = [sin(theta), -cos(theta)];

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% run_max = 5;

k = [10, 20, 40, 80, 160, 320, 640, 1280] %, 2560, 5120, 10240];

N_bf = [10, 20, 40, 80, 160, 320, 640, 1280];

filename{1} = 'poly_ss_k10_Nbf_10_1280.mat';
filename{2} = 'poly_ss_k20_Nbf_10_1280.mat';
filename{3} ='poly_ss_k40_Nbf_10_1280.mat';
filename{4} ='poly_ss_k80_Nbf_10_1280.mat';
filename{5} ='poly_ss_k160_Nbf_10_1280.mat';
filename{6} ='poly_ss_k320_Nbf_10_1280.mat';
filename{7} ='poly_ss_k640_Nbf_10_1280.mat';
filename{8} ='poly_ss_k1280_Nbf_10_1280.mat';

% 



for kk = 2:length(k)
    for n = 1:length(N_bf)
        tic
        % quadrature nodes and other information needed
        [G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, ...
            G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N,...
            G1_data.L] = discretistion_vars_graded(G1_data.G, C_wl, ...
            k(kk), Lgrad_coeff, alpha);
        
       
        [x1_1, y1_1, x1_2, y1_2, t1_bf_grid, t1_mid, ~, ~] = ...
        given_N_discretistion_vars_graded(G1_data.G, N_bf(n), Lgrad_coeff, alpha);
        
        % collocation points
        G1_data.s = [ t1_mid ; flip(G1_data.L - t1_mid) ];
        G1_data.x_col = [ x1_1 ; flip(x1_2) ];
        G1_data.y_col = [ y1_1 ; flip(y1_2) ];
        
            
        while length(G1_data.t_grid) < length(t1_bf_grid)  % case where there are more basis function than integration points, rendering the new bf useless
            C_wl= C_wl/2;
            % quadrature nodes and other information needed
            [G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N, G1_data.L] = ...
                discretistion_vars_graded(G1_data.G, C_wl, k(kk), Lgrad_coeff, alpha);
            
    %         [G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, G2_data.t_grid, G2_data.t_mid, G2_data.w, G2_data.N, G2_data.L] = ...
    %             discretistion_vars_graded(G2_data.G, C_wl, k(kk), Lgrad_coeff, alpha);
            disp('Increasing number of quadrature points due to basis functions')
        end
    
        [S11, uinc] = compute_matrix_single_screen_solve_graded(...
        G1_data, k(kk), t1_bf_grid, theta, C1, C2);
    
      
        
        coeffs = S11\uinc;
    
        t1_bf_grid_store{kk, n, :} = t1_bf_grid;
        
        
        aj_1{kk, n, :} = coeffs;

        C_wl_err = 1/40;
        
        [x1_1_q_err, y1_1_q_err, x1_2_q_err, y1_2_q_err, t1_grid_err, ...
            t1_mid_err, w1_err, N1_err, L1] = ...
            discretistion_vars_graded(G1_data.G, C_wl_err, k(kk), Lgrad_coeff, alpha);
        
        t1_mid_err = t1_mid_err(100:end);
        
        x1_plotting = [t1_mid_err; flip(L1 - t1_mid_err) ];
        
        w1_err = [w1_err(100:end); flip(w1_err(100:end))];
        
    %     phi_1 = zeros(length(N_bf), length(t1_mid_err)*2);
        
        phi_1{kk, n, :} = graded_coeff_2_solution(aj_1{kk, n}, t1_bf_grid_store{kk, n},...
                t1_mid_err, G1_data.L);
    
        G_mat = [G1_data.G(1) G1_data.G(2); G1_data.G(3) G1_data.G(4)];
    
        duidn_vec{kk, n, :} = 2*duidn(G_mat, G1_data.L, k(kk), d, x1_plotting);
    
%         figure()
%         plot(x1_plotting/G1_data.L, phi_1{kk, n, :}, 'DisplayName', ...
%             strcat('N = ', num2str(length(aj_1{kk, n}))))
%     
%         hold on
%         plot(x1_plotting/G1_data.L, duidn_vec{kk, n}, 'DisplayName', 'dui/dn')
%         xlabel('$x/L_1$')
%         ylabel('$\phi_1(x)$')
%         title(['Polynoial approximation of $\phi_{1}(x)$, when k = ', num2str(k(kk))])
%         legend show
%         xlim([-0.05 1.05])
%         ylim([-50 80])
%         fontsize(gca,18,"pixels")
%     
        phi_1_only_norm = sum(w1_err.*abs(duidn_vec{kk, n}));
    
        err_phi_1(kk, n) = sum(w1_err.*abs(phi_1{kk, n} - duidn_vec{kk, n}))./phi_1_only_norm;

        
        toc
        
    end
    warning('The saving of some information is not quite right')
    filename_str = filename{kk};

    N_bf_save = N_bf(n);

    aj_1_save = aj_1{kk, :};

    t1_bf_grid_save = t1_bf_grid_store{kk, :};

    err_phi_1_save = err_phi_1(kk, :);
% 
    save(filename_str, "G1_data", "N_bf_save", "aj_1_save", "t1_bf_grid_save", "err_phi_1_save")
% % % 
end

%% plotting

for kk = 1:length(k)

    for n = 1:length(N_bf)


%         C_wl_err = 1/40;
%         
%         [x1_1_q_err, y1_1_q_err, x1_2_q_err, y1_2_q_err, t1_grid_err, ...
%             t1_mid_err, w1_err, N1_err, L1] = ...
%             discretistion_vars_graded(G1_data.G, C_wl_err, k(kk), Lgrad_coeff, alpha);
%         
%         t1_mid_err = t1_mid_err(100:end);
%         
%         x1_plotting = [t1_mid_err; flip(L1 - t1_mid_err) ];
%         
%         w1_err = [w1_err(100:end); flip(w1_err(100:end))];
%         
%     %     phi_1 = zeros(length(N_bf), length(t1_mid_err)*2);
%         
%         phi_1{kk, n, :} = graded_coeff_2_solution(aj_1{kk, n}, t1_bf_grid_store{kk, n},...
%                 t1_mid_err, G1_data.L);
%     
%         G_mat = [G1_data.G(1) G1_data.G(2); G1_data.G(3) G1_data.G(4)];
%     
%         duidn_vec{kk, n, :} = 2*duidn(G_mat, G1_data.L, k(kk), d, x1_plotting);
    
        figure()
        plot(x1_plotting/G1_data.L, phi_1{kk, n, :}, 'DisplayName', ...
            strcat('N = ', num2str(length(aj_1{kk, n}))))
    
        hold on
        plot(x1_plotting/G1_data.L, duidn_vec{kk, n}, 'DisplayName', 'dui/dn')
        xlabel('$x/L_1$')
        ylabel('$\phi_1(x)$')
        title(['Polynoial approximation of $\phi_{1}(x)$, when k = ', num2str(k(kk))])
        legend show
        xlim([-0.05 1.05])
        ylim([-50 80])
        fontsize(gca,18,"pixels")
    
%         phi_1_only_norm = sum(w1_err.*abs(duidn_vec{kk, n}));
%     
%         err_phi_1(kk, n) = sum(w1_err.*abs(phi_1{kk, n} - duidn_vec{kk, n}))./phi_1_only_norm;

    end
        

end

err_phi_1





