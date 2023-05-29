% comparison between mine and Simons code

clear all

% adding paths
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/') % code to be tested
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Simon_code/') % simons code, taken to be as working code

%loading data, not sure if relevant but definatly will be helpful for pre
%determined information
load('newSCphi23_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')

theta = 0; %OC angle between the downwards verticle and incident wave, anti clockwise
d_test = [sin(theta), - cos(theta)];  % simon code, direction of incident wave

if d_test == d
    disp('direction of incident wave matches')
else
    error('d loaded does not match theta for thos case')
end
C_wl_true = C_wl(end);
clear C_wl
C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160];  % length of each interval C_wl a wavelength long

%%
%%% specific changes to make sure that code works regalrdless which way
%%% round the inputs of the vectors are, the comparing to high dof of SC
%%% won't now work
r1_startold = r1start;
r1_endold = r1end;
r2_startold = r2start;
r2_endold = r2end;
clear r1start r1end r2start r2end

% %case 1, switch G1 coordinates
% r1start = r1_endold;
% r1end = r1_startold;
% r2start = r2_startold;
% r2end = r2_endold;

% case 2, switch G2 coordinates
% r1start = r1_startold;
% r1end = r1_endold;
% r2start = r2_endold;
% r2end = r2_startold;
% 
% % case 3, switch G1 and G2 coordinates
r1start = r1_endold;
r1end = r1_startold;
r2start = r2_endold;
r2end = r2_startold;

%%
% specifics for my code
G1 = [r1start, r1end ];
G2 = [r2start, r2end ];

L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

C1 = 1;
C2 = pi;

% isolating the true high degrees of freedom solutions
phitrue1_0 = phi10;
phitrue2_1 = phi21;
phitrue1_2 = phi12;
phitrue2_3 = phi23;

% discretisation :
N1 = ceil(k*L1./(C_wl_true*2*pi)) % number of itervals on G1
N2 = ceil(k*L2./(C_wl_true*2*pi)) % number of intervals on G2

% simons solve
h1 = norm(r1end-r1start)/(N1+1); % the length of each element on screen 1
h2 = norm(r2end-r2start)/(N2+1); % the length of each element on screen 2
%     Next calculate the x and y coordinates of the element midpoints on screen
%     1 and then on screen 2

% x1_true = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
% y1_true = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;


x1_true = [0+h1:h1:L1-h1];
x2_true = [0+h2:h2:L2-h2];
for j = 1:length(C_wl)
    
    % discretisation :
    N1 = ceil(k*L1./(C_wl(j)*2*pi)) % number of itervals on G1
    N2 = ceil(k*L2./(C_wl(j)*2*pi)) % number of intervals on G2

    % simons solve
    h1 = norm(r1end-r1start)/N1; % the length of each element on screen 1
    h2 = norm(r2end-r2start)/N2; % the length of each element on screen 2
%     Next calculate the x and y coordinates of the element midpoints on screen
%     1 and then on screen 2
    
    x1 = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
    y1 = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;
    x2 = r2start(1)+((1:N2)-0.5)*(r2end(1)-r2start(1))/N2;
    y2 = r2start(2)+((1:N2)-0.5)*(r2end(2)-r2start(2))/N2;
    h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1
    h2vector = h2*ones(size(x2)); % and on screen 2
    h = [h1vector,h2vector]; % the lengths of the elements
    x = [x1,x2]; % the x-coords of the element midpoints
    y = [y1,y2]; % and their y-coords
    
    %step 0:
    phi1_0 = bem2(x1,y1,h1vector,k,d);
    
    %step 1:
    phi2_1 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_0,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
    
    %step 2:
    phi1_2 = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2_1,k,d);
    
    %step 3:
    phi2_3 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1_2,k,d);

    % oliver solve
    % method being tested OC:
    Q = 1;
    [A1, x_col1, A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
    [A2, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);
    
    [S21, col_points_1] = S21_op(G1, G2, k, N1, N2);
    [S12, col_points_2] = S21_op(G2, G1, k, N2, N1);
    
    
    % step 0
    u_i_1_0 = PW_incident(k, theta, G1, x_col1);
    aj_1_0 = A1\u_i_1_0.';
    
     % step 1
    u_i_2_1 =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_0, L1, N1-1, N1);
    aj_2_1 = A2\u_i_2_1;
    
    % step 2
    u_i_1_2 =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_1, L2, N2-1, N2);
    aj_1_2 = A1\u_i_1_2;
    
    % step 3
    u_i_2_3=  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_2, L1, N1-1, N1);
    aj_2_3= A2\u_i_2_3;
    
    
     % Now computing errors between the two, SC true same grid
    err_norm_SC_true_samegrid(j, 1) = norm(phi1_0 - aj_1_0, 1)/norm(phi1_0, 1);
    err_norm_SC_true_samegrid(j, 2) = norm(phi2_1 - aj_2_1, 1)/norm(phi2_1, 1);
    err_norm_SC_true_samegrid(j, 3) = norm(phi1_2 - aj_1_2, 1)/norm(phi1_2, 1);
    err_norm_SC_true_samegrid(j, 4) = norm(phi2_3 - aj_2_3, 1)/norm(phi2_3, 1);
    

    err_norm_eachhalftrue(j, 1) = norm(phi1_0 - aj_1_0, 1)/(norm(phi1_0, 1)/2+ norm(aj_1_0, 1)/2 );
    err_norm_eachhalftrue(j, 2) = norm(phi2_1 - aj_2_1, 1)/(norm(phi2_1, 1)/2+ norm(aj_2_1, 1)/2 );
    err_norm_eachhalftrue(j, 3) = norm(phi1_2 - aj_1_2, 1)/(norm(phi1_2, 1)/2+ norm(aj_1_2, 1)/2 );
    err_norm_eachhalftrue(j, 4) = norm(phi2_3 - aj_2_3, 1)/(norm(phi2_3, 1)/2+ norm(aj_2_3, 1)/2 );


    % fine grid mapped to coarse grid
    
    phi1_0true = interp1(x1_true, phitrue1_0, x_col1);
    phi2_1true = interp1(x2_true, phitrue2_1, x_col2);
    phi1_2true = interp1(x1_true, phitrue1_2, x_col1);
    phi2_3true = interp1(x2_true, phitrue2_3, x_col2);
    err_norm_SC_true_diffgrid(j, 1) = norm(phi1_0true - aj_1_0.', 1)/norm(phi1_0true, 1);
    err_norm_SC_true_diffgrid(j, 2) = norm(phi2_1true - aj_2_1.', 1)/norm(phi2_1true, 1);
    err_norm_SC_true_diffgrid(j, 3) = norm(phi1_2true - aj_1_2.', 1)/norm(phi1_2true, 1);
    err_norm_SC_true_diffgrid(j, 4) = norm(phi2_3true - aj_2_3.', 1)/norm(phi2_3true, 1);
    
    %plotting
    figure(1)
    plot(x_col1/L1, real(aj_1_0), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col1/L1, real(phi1_0), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    figure(2)
    plot(x_col2/L2, real(aj_2_1), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col2/L2, real(phi2_1), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    
    figure(3)
    plot(x_col1/L1, real(aj_1_2), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col1/L1, real(phi1_2), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    figure(4)
    plot(x_col2/L2, real(aj_2_3), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col2/L2, real(phi2_3), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))



end

figure(1)
title('Approximations for $\phi_{1}^{(0)}$')
legend show

figure(2)
title('Approximations for $\phi_{2}^{(1)}$')
legend show

figure(3)
title('Approximations for $\phi_{1}^{(2)}$')
legend show

figure(4)
title('Approximations for $\phi_{2}^{(3)}$')
legend show

%%
r_plotting = [0:1:3];
%plotting
figure()
hold on
for j = 1:length(C_wl)
    
    semilogy(r_plotting, err_norm_SC_true_samegrid(j, :), 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )

end


%%
% EOC
for j = 1:length(C_wl) - 1
    for r = 1:4
    
    EOC_samegrid_SCtrue(j, r) = log2(err_norm_SC_true_samegrid(j, r)/err_norm_SC_true_samegrid(j+1, r));
    
    EOC_samegrid_bothtrue(j, r) = log2(err_norm_eachhalftrue(j, r)/err_norm_eachhalftrue(j+1, r));
    
    EOC_diffgrid_SCtrue(j, r) = log2(err_norm_SC_true_diffgrid(j, r)/err_norm_SC_true_diffgrid(j+1, r));
    
    end
end
