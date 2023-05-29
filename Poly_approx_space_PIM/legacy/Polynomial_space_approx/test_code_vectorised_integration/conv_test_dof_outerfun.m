% this script will compute the scattered field from two disjoint screens in
% any orientation in a polynomial appproaximation space. The singularities
% are handled by using a product integration method and all integrations
% are approximated using the midpoint rule.
%
% In this script we will be testing whether we have convergence for
% various different point in each of the matrices computed.

clear all

% adding path for Simon's code
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Simon_code/') % simons code, taken to be as working code

tic

% parallel screens
% len = 4*pi; % length ish of screen
% D = 2*pi + 0.005; % distance between screens

% G1 = [-len, len, 0, 0];
% G2 = [-len+D, len, D , 0];
% introducing the screens
% Introducing the screens
G1 = [-2*pi, 2*pi, 0, 0];
% G1 = [-pi, 2*pi, 0, 0];
% G1 = [0, 0, -2*pi, 2*pi];
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1

G2 = [ pi, 0, 3*pi, 2*pi];
% G2 = [pi, 0, 2*pi, 3*pi];
% G2 = [5*pi, 3*pi, 2*pi, 0];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

% C_wl= [1/5, 1/10, 1/20, 1/40, 1/80, 1/160, 1/320];

C_wl = 1/20;

k = 10;  % wavenumber

R_approx = 2;

f = @(x) 1;

theta = 0;

C1 = 1; C2 = pi;

% adapting my inputs to Simons
r1start = [G1(1), G1(2)]; r1end = [G1(3), G1(4)];
r2start = [G2(1), G2(2)]; r2end = [G2(3), G2(4)];
d = [sin(theta), -cos(theta)];

% computing the "true" solution, large number of degrees of freedom, large
% reflection order
% C_wl_true = C_wl(end)*2;
R_true = 3*R_approx;
% [x1, y1, t1, t1_mid, h1, hvector1, N1, L1] = discretisation_variables(G1, C_wl, k);
% [x2, y2, t2, t2_mid, h2, hvector2, N2, L2] = discretisation_variables(G2, C_wl, k);
% [aj_1_r_true, aj_2_r_true, ~, ~, ~, ~] = iterative_2screen_solver(k, G1, G2, f, theta, C_wl, R_true, C1, C2);
%     

% looping over the number of degrees of freedom and also convergence with
% respect to iterations
for j = 1:length(C_wl)
    [x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, C_wl(j), k);
    [x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, C_wl(j), k);
    
    % solving for r reflection
    [aj_1_r, aj_2_r, A1, A2, S12, S21] = iterative_2screen_solver(k, G1, G2, f, theta, C_wl(j), R_true, C1, C2);
%     
%     aj_1_test = aj_1_r(:, :, j);
%     aj_2_test = aj_2_r(:, :, j);
    
    % constructing and doing the all in one solve
    A = [A1 S12;
        S21 A2];
    u_i = [robhop_PW_incident(k, theta, x1, y1) robhop_PW_incident(k, theta, x2, y2)];
    
    aj_allin1_1_2 = A\u_i.';
    
    aj_1_allone = aj_allin1_1_2(1:length(t1_mid), 1);
    aj_2_allone = aj_allin1_1_2(length(t1_mid)+1: length(t1_mid) + length(t2_mid), 1);
    
    % simons code to compare to
%     [phi1, b] = bem2(x1,y1,h1vector,k,d);
%     SC_OP_ui_comp = [b robhop_PW_incident(k, theta, x1, y1).']
%     
%     phi10 = phi1; %Keeping the zero iteration
% 
%     phi2 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
%     phi = [phi1.',phi2.'];
%     phi21 = phi2; % keeping the 1 iteration
%     
%     phi1 = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phi2,k,d);
%     phi = [phi1.',phi2.'];
%     phi12 = phi1; 
%     
%     phi2 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phi1,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
%     phi = [phi1.',phi2.'];
%     phi23 = phi2;
    
    
        
%     err_SC_true(j, 1) = norm( phi10 - aj_1_r(:, 1))/norm( phi10);
%     err_SC_true(j, 2) = norm( phi21 - aj_2_r(:, 1))/norm( phi21);
%     err_SC_true(j, 3) = norm( phi12 - aj_1_r(:, 2))/norm( phi12);
%     err_SC_true(j, 4) = norm( phi23 - aj_2_r(:, 2))/norm( phi23);   

    for r = 1:R_approx
        err_G1_wrtR(r, j) = norm(aj_1_r(:, end) - aj_1_r(:, r), 1 )/norm(aj_1_r(:, end), 1);
        
        err_G2_wrtR(r, j) = norm(aj_2_r(:, end) - aj_2_r(:, r) , 1)/norm(aj_2_r(:, end), 1);
        
        err_G1_all1_true(r, j) = norm( aj_1_allone - aj_1_r(:, r) , 1)/norm(aj_1_allone, 1);
        err_G2_all1_true(r, j) = norm( aj_2_allone - aj_2_r(:, r), 1 )/norm(aj_2_allone, 1);
        
        if r == 1
            
        else
            err_G1_previous(r, j) = norm(aj_1_r(:, r) - aj_1_r(:, r-1), 1)/norm(aj_1_r(:, r), 1);
            err_G2_previous(r, j) = norm(aj_2_r(:, r) - aj_2_r(:, r-1), 1)/norm(aj_2_r(:, r), 1);
        end 
        
    end
end

%% plotting errors
% convergence wrt number of iterations
r1 = [0:2:R_approx*2 - 2];
r2 = [1:2:R_approx*2 - 1];
figure()
for j = 1:length(C_wl)
    semilogy(r1, err_G1_wrtR(:, j), 'DisplayName', sprintf('dof per wl = %g', 1/C_wl(j)) )
    hold on
end
xlabel('r')
ylabel('$\Vert \phi_{1}^{R} - \phi_{1}^{r} \Vert / \Vert \phi_{1}^{R} \Vert$')
title('Error plot on $\Gamma_{1}$ with respect to number of iterations compared to high order of iterations')
legend show

figure()
for j = 1:length(C_wl)
    semilogy(r2, err_G2_wrtR(:, j), 'DisplayName', sprintf('dof per wl = %g', 1/C_wl(j)) )
    hold on
end
xlabel('r')
ylabel('$\Vert \phi_{2}^{R} - \phi_{2}^{r} \Vert / \Vert \phi_{2}^{R} \Vert$')
title('Error plot on $\Gamma_{2}$ with respect to number of iterations compared to high order of iterations')
legend show

% plotting errors where one step solution is taken as true
figure()
for j = 1:length(C_wl)
    semilogy(r1, err_G1_all1_true(:, j), 'DisplayName', sprintf('dof per wl = %g', 1/C_wl(j)) )
    hold on
end
xlabel('r')
ylabel('$\Vert \tilde{\phi_{1}} - \phi_{1}^{r} \Vert / \Vert \tilde{\phi_{1}} \Vert$')
title('Error plot on $\Gamma_{1}$ with respect to number of iterations compared to a standard BEM direct solver')
legend show

figure()
for j = 1:length(C_wl)
    semilogy(r2, err_G2_all1_true(:, j), 'DisplayName', sprintf('dof per wl = %g', 1/C_wl(j)) )
    hold on
end
xlabel('r')
ylabel('$\Vert \tilde{\phi_{2}} - \phi_{2}^{r} \Vert / \Vert \tilde{\phi_{2}} \Vert$')
title('Error plot on $\Gamma_{2}$ with respect to number of iterations compared to a standard BEM direct solver')
legend show


%%
% plotting for specific # dof per wavelength
for r = 1:R_approx
    
    figure(50)
    plot(t1_mid/L1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
    hold on
    
    figure(51)
    plot(t2_mid/L2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 1))
    hold on

end
figure(50)
% plot(t1_mid/L1, real(phi12), 'DisplayName', 'SC r = 2')
title('Solution on $\Gamma_{1}$ for an increasing number of iterations')
xlim([ -0.05 1.05 ])
ylim([-20 80])
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{r}(x)$')
legend show

figure(51)
% plot(t2_mid/L2, real(phi23), 'DisplayName', 'SC r = 3')
title('Solution on $\Gamma_{2}$ for an increasing number of iterations')
xlim([ -0.05 1.05 ])
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{r}(x)$')
legend show

%% plotting in the domain
keyboard
% now that we have the coefficients we can plot in the domain
% Makr a square grid around the screens
X_coordinates = [G1(1), G1(3), G2(1), G2(3)];
X_min = min(X_coordinates);
X_max = max(X_coordinates);
X_diff = abs(X_max - X_min);

Y_coordinates = [G1(2), G1(4), G2(2), G2(4)];
Y_min = min(Y_coordinates);
Y_max = max(Y_coordinates);
Y_diff = abs(Y_max - Y_min);

% maybe also change this part
N_diff = max(abs(Y_diff - X_diff)) + 5*pi;  % largest distance either x direction or y + 4*pi (2*pi in either direction)

N_d = 500;

min_leftorbottom = min(X_min, Y_min);
max_toporright = max(X_max, Y_max);
h_d = N_diff/N_d;
% need to change this bit so that it is closer zoomed in on the picture
X = [min_leftorbottom - 5: h_d: max_toporright + 5];
Y = [min_leftorbottom - 5: h_d: max_toporright + 5];
% first compute the incident field


% incident wave:
U_i = incident_domain(k, theta, X, Y);

figure(); pcolor(X, Y, real(U_i)); shading interp

% step r = 0:
tic
[us1_0, x1, y1] = scattered_domain( G1, aj_1_r(:, 1), k, X, Y, N1);
figure(); pcolor(X, Y, real(us1_0)); shading interp; colorbar; title('scattered field no-rereflections ')
toc

UT_1_0 = U_i - us1_0;

figure()
pcolor(X, Y, real(UT_1_0));
shading interp;
colormap(jet);
title('Total field after 0 reflection')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];


tic
[us2_1, x1, y1, x2, y2] = scattered_domain2( G1, G2, aj_1_r(:, 1), aj_2_r(:, 1), k, X, Y, N1, N2);
toc

UT_2_1 = U_i - us2_1;

figure(); pcolor(X, Y, real(UT_2_1)); shading interp; colormap(jet); colorbar; title('total field 1 rereflections ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

figure()
    pcolor(X, Y, real(UT_2_1 - UT_1_0));
    shading interp;
    colormap(jet);
    title(['Difference between total field of ',num2str(1),' reflections and ',num2str(0),' reflections'])
    hold on
    Gamma_1 = plot(x1, y1);
    Gamma_1.LineWidth = 4;
    Gamma_1.Color = [0 0 0];
    Gamma_2 = plot(x2, y2);
    Gamma_2.LineWidth = 4;
    Gamma_2.Color = [0 0 0];

UT_2_r = UT_2_1;

for j = 2:4
    
    %solution on \Gamma_{1}
    [us1_r, ~, ~] = scattered_domain2( G1, G2, aj_1_r(:, j), aj_2_r(:, j-1), k, X, Y, N1, N2);
    UT_1_r = U_i - us1_r;
    
    figure()
    pcolor(X, Y, real(UT_1_r));
    shading interp;
    colormap(jet);
    title(['Total field after ',num2str(2*j - 2),' reflections'])
    hold on
    Gamma_1 = plot(x1, y1);
    Gamma_1.LineWidth = 4;
    Gamma_1.Color = [0 0 0];
    Gamma_2 = plot(x2, y2);
    Gamma_2.LineWidth = 4;
    Gamma_2.Color = [0 0 0];
    
    figure()
    pcolor(X, Y, real(UT_1_r - UT_2_r));
    shading interp;
    colormap(jet);
    title(['Difference between total field of ',num2str(2*j -2),' reflections and ',num2str(2*j-3),' reflections'])
    hold on
    Gamma_1 = plot(x1, y1);
    Gamma_1.LineWidth = 4;
    Gamma_1.Color = [0 0 0];
    Gamma_2 = plot(x2, y2);
    Gamma_2.LineWidth = 4;
    Gamma_2.Color = [0 0 0];
    
    %solution on \Gamma_{2}
    [us2_r, ~, ~] = scattered_domain2( G1, G2, aj_1_r(:, j), aj_2_r(:, j), k, X, Y, N1, N2);
    UT_2_r = U_i - us2_r;
    
    figure()
    pcolor(X, Y, real(UT_2_r));
    shading interp;
    colormap(jet);
    title(['Total field after ',num2str(2*j - 1),' reflections'])
    hold on
    Gamma_1 = plot(x1, y1);
    Gamma_1.LineWidth = 4;
    Gamma_1.Color = [0 0 0];
    Gamma_2 = plot(x2, y2);
    Gamma_2.LineWidth = 4;
    Gamma_2.Color = [0 0 0];
    
    figure()
    pcolor(X, Y, real(UT_2_r - UT_1_r));
    shading interp;
    colormap(jet);
    title(['Difference between total field of ',num2str(2*j - 1),' reflections and ',num2str(2*j-2),' reflections'])
    hold on
    Gamma_1 = plot(x1, y1);
    Gamma_1.LineWidth = 4;
    Gamma_1.Color = [0 0 0];
    Gamma_2 = plot(x2, y2);
    Gamma_2.LineWidth = 4;
    Gamma_2.Color = [0 0 0];

    

    
end

keyboard
%% plotting high order plot in domain
% now that we have the coefficients we can plot in the domain
% Makr a square grid around the screens
X_coordinates = [G1(1), G1(3), G2(1), G2(3)];
X_min = min(X_coordinates);
X_max = max(X_coordinates);
X_diff = abs(X_max - X_min);

Y_coordinates = [G1(2), G1(4), G2(2), G2(4)];
Y_min = min(Y_coordinates);
Y_max = max(Y_coordinates);
Y_diff = abs(Y_max - Y_min);

% maybe also change this part
N_diff = max(abs(Y_diff - X_diff)) + 5*pi;  % largest distance either x direction or y + 4*pi (2*pi in either direction)

N_d = 500;

min_leftorbottom = min(X_min, Y_min);
max_toporright = max(X_max, Y_max);
h_d = N_diff/N_d;
% need to change this bit so that it is closer zoomed in on the picture
X = [min_leftorbottom - 5: h_d: max_toporright + 5];
Y = [min_leftorbottom - 5: h_d: max_toporright + 5];

% incident wave:
U_i = incident_domain(k, theta, X, Y);

% scattered field
tic
[us2_highorder, x1, y1, x2, y2] = scattered_domain2( G1, G2, aj_1_r(:, end), aj_2_r(:, end), k, X, Y, N1, N2);
toc

U_T_highorder = U_i - us2_highorder;
figure(); pcolor(X, Y, real(U_T_highorder)); shading interp; colormap(jet); colorbar; 
% title('total field 1 rereflections ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

