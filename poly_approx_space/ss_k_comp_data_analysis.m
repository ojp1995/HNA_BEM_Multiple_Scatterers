% single screen k convergence reading in data and then producing plots etc

clear all

addpath('../General_functions/')

% begin reading in functions
filename{1} = 'poly_ss_k10_Nbf_10_1280.mat';
filename{2} = 'poly_ss_k20_Nbf_10_1280.mat';
filename{3} ='poly_ss_k40_Nbf_10_1280.mat';
filename{4} ='poly_ss_k80_Nbf_10_1280.mat';
filename{5} ='poly_ss_k160_Nbf_10_1280.mat';
filename{6} ='poly_ss_k320_Nbf_10_1280.mat';
filename{7} ='poly_ss_k640_Nbf_10_1280.mat';
filename{8} ='poly_ss_k1280_Nbf_10_1280.mat';

% keyboard

for j = 1:length(filename)
    
    tmp_file = open(filename{j});

    err_phi_1(j, :) = tmp_file.err_phi_1_save;

    clear tmp_file

end

k = [10, 20, 40, 80, 160, 320, 640, 1280];
N_bf = [10, 20, 40, 80, 160, 320, 640, 1280];

figure(); 
for j = 1:length(k)
    loglog(N_bf(1:end), err_phi_1(j, :), 'DisplayName', strcat('k = ', num2str(k(j))))
    hold on
end

legend show
xlabel('Number of degrees of freedom')
ylabel('Relative $L^{1}$ error ')
title('Relative $L^{1}$ error comparing our approximation to the analytic solution as $k \rightarrow \infty$.')
fontsize(gca,18,"pixels")