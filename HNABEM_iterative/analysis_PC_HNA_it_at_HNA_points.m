% comparison between the PC and HNA iterative method evaluating at the
% points that the HNA solution was approximated at.

clear all

% adding paths

addpath('../General_functions/')
addpath('../../BEAM_HNABEMLAB/')
addPathsHNA  % allows HNABEM to find all of the relevatn subfolders

% addpath('test5a_HF_data_p4_8/')
% loading in data from HNA solver

test = {'test1_HNA_pmax4_overlap2',...
    'test1_HNA_pmax5_overlap2',...
    'test1_HNA_pmax6_overlap2', ...
    'test1_HNA_pmax7_overlap2', ...
    'test1_HNA_pmax8_overlap2'}
% test = {'test2_HNA_pmax4_overlap2_pi_4',...
%     'test2_HNA_pmax5_overlap2_pi_4',...
%     'test2_HNA_pmax6_overlap2_pi_4', ...
%     'test2_HNA_pmax7_overlap2_pi_4', ...
%     'test2_HNA_pmax8_overlap2_pi_4'}
% test = {'test3_HNA_pmax4_overlap2',...
%     'test3_HNA_pmax5_overlap2',...
%     'test3_HNA_pmax6_overlap2', ...
%     'test3_HNA_pmax7_overlap2', ...
%     'test3_HNA_pmax8_overlap2'}
% test = {'test4_HNA_pmax4_overlap2',...
%     'test4_HNA_pmax5_overlap2',...
%     'test4_HNA_pmax6_overlap2', ...
%     'test4_HNA_pmax7_overlap2', ...
%     'test4_HNA_pmax8_overlap2'}

% addpath('test5a_HF_data_p4_8/')
% test = {'test5a_HNA_pmax4_overlap2_dof5_80',...
%     'test5a_HNA_pmax5_overlap2_dof5_80.mat', ...
%     'test5a_HNA_pmax6_overlap2_dof5_80.mat', ...
%     'test5a_HNA_pmax7_overlap2_dof5_80.mat', ...
%     'test5a_HNA_pmax8_overlap2_dof5_80.mat'}


% test = {'test5b_HNA_pmax4_overlap2_pi_4_short.mat',...
%     'test5b_HNA_pmax5_overlap2_pi_4_short.mat', ...
%     'test5b_HNA_pmax6_overlap2_pi_4_short.mat', ...
%     'test5b_HNA_pmax7_overlap2_pi_4_short.mat', ...
%     'test5b_HNA_pmax8_overlap2_pi_4_short.mat'}


%% loading in data from piecewise constant direct solver
 
addpath('../poly_approx_space/')
load('PC_direct_test1_k10_theta0.mat')
% load('PC_direct_test2_k10_theta_pi_4.mat')
% load('PC_direct_test3_k10_theta0.mat')
% load('PC_direct_test4_k10_theta0.mat')
% load('PC_direct_test5a_k10_thetapi_4.mat')
% load('PC_direct_test5b_k10_thetapi_4_short.mat')
%%
G1_data_poly.G = info_needed.G1;
G2_data_poly.G = info_needed.G2;



k = info_needed.k;
theta = info_needed.theta;
Lgrad_coeff_poly = info_needed.L_grad_coeff;
alpha_poly = info_needed.alpha;
bf_dof_per_wl = info_needed.bf_dof_per_wl;

for l = 1:length(test) % looping through tests

    load(test{l})

    for j = 1:length(phi1_HNA)
    
        for n = 1:length(aj1_coeff)
    
            G1_data_poly =  get_bf_graded_grid(G1_data_poly, bf_dof_per_wl(n), k, ...
                    Lgrad_coeff_poly, alpha_poly);
                
            G2_data_poly =  get_bf_graded_grid(G2_data_poly, bf_dof_per_wl(n), k, ...
                    Lgrad_coeff_poly, alpha_poly);
        
            phi_1_poly = graded_coeff_2_solution(aj1_coeff{n}, ...
                    G1_data_poly.t_bf_grid, G1_data_HNA{j}.t_mid_q_outer, ...
                    G1_data_poly.L);
        
            phi_2_poly = graded_coeff_2_solution(aj2_coeff{n}, ...
                    G2_data_poly.t_bf_grid, G2_data_HNA{j}.t_mid_q_outer, ...
                    G2_data_poly.L);
        
            err_1(n, j) = sum((abs(phi_1_poly - phi1_HNA{j}{end})./...
                abs(phi_1_poly)).*G1_data_HNA{j}.w_comb_outer);

            err_1_test(n, j) = sum(abs(phi_1_poly - phi1_HNA{j}{end}).*G1_data_HNA{j}.w_comb_outer)/...
                sum(abs(phi_1_poly).*G1_data_HNA{j}.w_comb_outer);
        
            err_2(n, j) = sum((abs(phi_2_poly - phi2_HNA{j}{end})./...
                    abs(phi_2_poly)).*G2_data_HNA{j}.w_comb_outer);

            err_2_test(n, j) = sum(abs(phi_2_poly - phi2_HNA{j}{end}).*G2_data_HNA{j}.w_comb_outer)/...
                sum(abs(phi_2_poly).*G2_data_HNA{j}.w_comb_outer);
    
    
        end
    end


    err1_diag(l, :) = diag(err_1).';

    err1_diag_test(l, :) = diag(err_1_test).';

    err2_diag(l, :) = diag(err_2).';

    err2_diag_test(l, :) = diag(err_2_test).';

%     phi1_dof(l) = length(phi1_HNA{l}{end});
%     phi2_dof(l) = length(phi2_HNA{l}{end});
end

% v_N1_HNA_cell{end}{end}.Ndim, err1_diag
% v_N2_HNA_cell{end}{end}.Ndim, err2_diag

% %%
% figure();
% for j = 1:length(bf_dof_per_wl)
%     txt = ['bf dof = ', mat2str(phi1_dof(j))];
%     semilogy(1./bf_dof_per_wl, err1_diag(j, :), 'DisplayName', txt)
%     hold on
% 
% end
% xlabel('Quadrature points per wavelength')
% ylabel('$L^{1}$ relative error')
% title('Comparison between piecewise constant conventional BEM and HNA iterative method approximation of $\phi_{1}$')
% legend show

%%
% keyboard
% computed manually for test 5a
phi1_dof = [64, 90, 118, 152, 188];
phi2_dof = [64, 90, 118, 152,188];
figure();
for j = 1:length(bf_dof_per_wl)
    txt = ['quad points per wavelength = ', mat2str(1/bf_dof_per_wl(j))];
    semilogy(phi1_dof, err1_diag(:, j), 'DisplayName', txt)
    hold on

end
xlabel('Degrees of freedom')
ylabel('$L^{1}$ relative error')
title('Comparison between piecewise constant conventional BEM and HNA iterative method approximation of $\phi_{1}$')
legend show
ylim([5e-3 2e0])

figure();
for j = 1:length(bf_dof_per_wl)
    txt = ['quad points per wavelength = ', mat2str(1/bf_dof_per_wl(j))];
    semilogy(phi1_dof, err1_diag_test(:, j), 'DisplayName', txt)
    hold on

end
xlabel('Degrees of freedom')
ylabel('$L^{1}$ relative error')
title('Comparison between piecewise constant conventional BEM and HNA iterative method approximation of $\phi_{1}$')
legend show
ylim([ 5e-4 1e0 ])



figure();
for j = 1:length(bf_dof_per_wl)
    txt = ['quad points per wavelength = ', mat2str(1/bf_dof_per_wl(j))];
    semilogy(phi2_dof, err2_diag(:, j), 'DisplayName', txt)
    hold on

end
xlabel('Degrees of freedom')
ylabel('$L^{1}$ relative error')
title('Comparison between piecewise constant conventional BEM and HNA iterative method approximation of $\phi_{2}$')
legend show
ylim([ 5e-4 1e-1 ])
ylim([1e-3 1e0])

figure();
for j = 1:length(bf_dof_per_wl)
    txt = ['quad points per wavelength = ', mat2str(1/bf_dof_per_wl(j))];
    semilogy(phi2_dof, err2_diag_test(:, j), 'DisplayName', txt)
    hold on

end
xlabel('Degrees of freedom')
ylabel('$L^{1}$ relative error')
title('Comparison between piecewise constant conventional BEM and HNA iterative method approximation of $\phi_{2}$')
legend show
ylim([ 5e-4 1e0 ])
% ylim([ 1e-9 1e-3  ])

