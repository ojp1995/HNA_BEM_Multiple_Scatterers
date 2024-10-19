% increasing k comparison between conventional BEM and Hf iterative method

clear all

addpath('../General_functions/')
% addpath('../../BEAM_HNABEMLAB/')

% add path to data HF iterative method
addpath('HF_it_k_increase/')
% loading data for each individual k from relevant folder
 k1 = {'HF_it_test1_k1_dof40_theta0_short_p4', ...
    'HF_it_test1_k1_dof40_theta0_short_p5', ...
    'HF_it_test1_k1_dof40_theta0_short_p6', ...
    'HF_it_test1_k1_dof40_theta0_short_p7', ...
    'HF_it_test1_k1_dof40_theta0_short_p8'};

k2 = {'HF_it_test1_k2_dof40_theta0_short_p4', ...
    'HF_it_test1_k2_dof40_theta0_short_p5', ...
    'HF_it_test1_k2_dof40_theta0_short_p6', ...
    'HF_it_test1_k2_dof40_theta0_short_p7', ...
    'HF_it_test1_k2_dof40_theta0_short_p8'};

k5 = {'HF_it_test1_k5_dof40_theta0_short_p4', ...
    'HF_it_test1_k5_dof40_theta0_short_p5', ...
    'HF_it_test1_k5_dof40_theta0_short_p6', ...
    'HF_it_test1_k5_dof40_theta0_short_p7', ...
    'HF_it_test1_k5_dof40_theta0_short_p8'};

k10 = {'HF_it_test1_k10_dof40_theta0_short_p4', ...
    'HF_it_test1_k10_dof40_theta0_short_p5', ...
    'HF_it_test1_k10_dof40_theta0_short_p6', ...
    'HF_it_test1_k10_dof40_theta0_short_p7', ...
    'HF_it_test1_k10_dof40_theta0_short_p8'};

k20 = {'HF_it_test1_k20_dof40_theta0_short_p4', ...
    'HF_it_test1_k20_dof40_theta0_short_p5', ...
    'HF_it_test1_k20_dof40_theta0_short_p6', ...
    'HF_it_test1_k20_dof40_theta0_short_p7', ...
    'HF_it_test1_k20_dof40_theta0_short_p8'};

k40 = {'HF_it_test1_k40_dof40_theta0_short_p4', ...
    'HF_it_test1_k40_dof40_theta0_short_p5', ...
    'HF_it_test1_k40_dof40_theta0_short_p6', ...
    'HF_it_test1_k40_dof40_theta0_short_p7', ...
    'HF_it_test1_k40_dof40_theta0_short_p8'};

HF_test = {k1, k2, k5, k10, k20, k40};

% loading PC data into a cell
addpath('../poly_approx_space/PC_k_increase_test/')
PC_test = {'PC_direct_test1_k1_dof40_theta0_short', ...
    'PC_direct_test1_k2_dof40_theta0_short', ...
    'PC_direct_test1_k5_dof40_theta0_short', ...
    'PC_direct_test1_k10_dof40_theta0_short', ...
    'PC_direct_test1_k20_dof40_theta0_short', ...
    'PC_direct_test1_k40_dof40_theta0_short'};


for k_test = 1:length(PC_test)  % looping through for each k experiment
    load(PC_test{k_test});
    
    % extracting import info 
    G1_data_poly.G = info_needed.G1;
    G2_data_poly.G = info_needed.G2;
 
    k = info_needed.k;
    theta = info_needed.theta;
    Lgrad_coeff_poly = info_needed.L_grad_coeff;
    alpha_poly = info_needed.alpha;
    bf_dof_per_wl = info_needed.bf_dof_per_wl;

    for hf = 1:length(k10)  % looping through each of the individual tests, i.e., p = 4, 5, 6, 7, 8
        load(HF_test{k_test}{hf});

        G1_data_poly =  get_bf_graded_grid(G1_data_poly, bf_dof_per_wl, k, ...
                    Lgrad_coeff_poly, alpha_poly);
                
        G2_data_poly =  get_bf_graded_grid(G2_data_poly, bf_dof_per_wl, k, ...
                    Lgrad_coeff_poly, alpha_poly);

        phi_1_poly = graded_coeff_2_solution(aj1_coeff{1}, ...
                    G1_data_poly.t_bf_grid, G1_data_HNA{1}.t_mid_q_outer, ...
                    G1_data_poly.L);
        
        phi_2_poly = graded_coeff_2_solution(aj2_coeff{1}, ...
                    G2_data_poly.t_bf_grid, G2_data_HNA{1}.t_mid_q_outer, ...
                    G2_data_poly.L);

        err_1(k_test, hf) = sum((abs(phi_1_poly - phi1_HNA{1}{end})./...
                abs(phi_1_poly)).*G1_data_HNA{1}.w_comb_outer);
        
        err_1_test(k_test, hf) = sum(abs(phi_1_poly - phi1_HNA{1}{end}).*G1_data_HNA{1}.w_comb_outer)/...
                sum(abs(phi_1_poly).*G1_data_HNA{1}.w_comb_outer); 

        err_2(k_test, hf) = sum((abs(phi_2_poly - phi2_HNA{1}{end})./...
                    abs(phi_2_poly)).*G2_data_HNA{1}.w_comb_outer);

        err_2_test(k_test, hf) = sum(abs(phi_2_poly - phi2_HNA{1}{end}).*G2_data_HNA{1}.w_comb_outer)/...
                    sum(abs(phi_2_poly).*G2_data_HNA{1}.w_comb_outer);

%         phi1_dof(k_test, hf) = v_N1_HNA_cell{1}{end}.Ndim;
%         phi2_dof(k_test, hf) = v_N2_HNA_cell{1}{end}.Ndim;
%     
    end 

    



end

Ndof = [64, 90, 118, 152, 188];
k_vals = [1, 2, 5, 10, 20, 40];

err_1, err_2 % phi1_dof, phi2_dof

%%
% plotting the error with each k being on a different line. Error as a
% function as dof
figure()
for k_plot = 1:length(PC_test)
    txt = ['k = ', mat2str(k_vals(k_plot))];
    semilogy(Ndof, err_1(k_plot, :), 'DisplayName', txt)
    hold on

end
xlabel('degrees of freedom of HNA method')
ylabel('Relative $L^{1}$ error')
title('Relative error for high-frequency iterative method for increasing $k$ of $\phi_{1}$')
legend show
ylim([1e-3 1e-1])

figure()
for k_plot = 1:length(PC_test)
    txt = ['k = ', mat2str(k_vals(k_plot))];
    semilogy(Ndof, err_1_test(k_plot, :), 'DisplayName', txt)
    hold on

end
xlabel('degrees of freedom of HNA method')
ylabel('Relative $L^{1}$ error')
title('Relative error for high-frequency iterative method for increasing $k$ of $\phi_{1}$')
legend show

figure()
for k_plot = 1:length(PC_test)
    txt = ['k = ', mat2str(k_vals(k_plot))];
    semilogy(Ndof, err_2(k_plot, :), 'DisplayName', txt)
    hold on

end
xlabel('degrees of freedom of HNA method')
ylabel('Relative $L^{1}$ error')
title('Relative error for high-frequency iterative method for increasing $k$ of $\phi_{2}$')
legend show
ylim([1e-3 1e-1])


figure()
for k_plot = 1:length(PC_test)
    txt = ['k = ', mat2str(k_vals(k_plot))];
    semilogy(Ndof, err_2_test(k_plot, :), 'DisplayName', txt)
    hold on

end
xlabel('degrees of freedom of HNA method')
ylabel('Relative $L^{1}$ error')
title('Relative error for high-frequency iterative method for increasing $k$ of $\phi_{2}$')
legend show

