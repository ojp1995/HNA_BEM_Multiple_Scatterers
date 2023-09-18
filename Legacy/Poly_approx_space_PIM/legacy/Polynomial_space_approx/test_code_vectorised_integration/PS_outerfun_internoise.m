% in this script we will be attempting to have not 1 but 2 point sources,
% idea is motivation for sound barriers

clear all

tic

% simon noise barrier paper
xps1 = 8;
yps1 = 0.5;

xps2 = 8;
yps2 = -0.5;
% screens
G1 = [0, -3, 0, 3];
G2 = [34, -3, 34, 3];

% H = 3*pi; % Height of each tip of screen from y = 0, i.e screens are 2H high with centre at y = 0;
% dist = 2*pi; % distance each screen is from the the point source
% 
% xps1 = 0;
% yps1 = H/2;
% 
% xps2 = 0;
% yps2 = -H/2;



% G1 = [-dist, -H, -dist, H];
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1


% G2 = [dist, -H, dist, H];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

C_wl= [1/5, 1/10, 1/20, 1/40] %, 1/80, 1/160] %, 1/320];

% C_wl = 1/20;

k = 2*pi/0.72886;  % wavenumber

R_approx = 60;

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

for j = 1:length(C_wl)
    [x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, C_wl(j), k);
    [x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, C_wl(j), k);
    
    % solving for r reflection
    [aj_1_r, aj_2_r, A1, A2, S12, S21] = PS_2_iterative_2screen_solver(k, G1, G2, f, xps1, yps1, xps2, yps2, C_wl(j), R_approx, C1, C2);

%     aj_1_test = aj_1_r(:, :, j);
%     aj_2_test = aj_2_r(:, :, j);
    
    % constructing and doing the all in one solve
    A = [A1 S12;
        S21 A2];
    u_i = [(PS_incident(k, xps1, yps1, x1, y1) + PS_incident(k, xps2, yps2, x1, y1)) PS_incident(k, xps1, yps1, x2, y2) + PS_incident(k, xps2, yps2, x2, y2)];
    
    aj_allin1_1_2 = A\u_i.';
    
    aj_1_allone = aj_allin1_1_2(1:length(t1_mid), 1);
    aj_2_allone = aj_allin1_1_2(length(t1_mid)+1: length(t1_mid) + length(t2_mid), 1);
    
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
for r = 1:10:R_approx
    
    figure(52)
    plot(t1_mid/L1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
    hold on
    
    figure(53)
    plot(t2_mid/L2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 1))
    hold on

end
figure(52)
% plot(t1_mid/L1, real(phi12), 'DisplayName', 'SC r = 2')
title('Solution on $\Gamma_{1}$ for an increasing number of iterations')
xlim([ -0.05 1.05 ])
% ylim([-20 80])
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{r}(x)$')
legend show

figure(53)
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
N_diff = max(abs(Y_diff - X_diff)) + 2*pi;  % largest distance either x direction or y + 4*pi (2*pi in either direction)

N_d = 300;

min_leftorbottom = min(X_min, Y_min);
max_toporright = max(X_max, Y_max);
h_d = N_diff/N_d;
% need to change this bit so that it is closer zoomed in on the picture
% X = [min_leftorbottom - 5: h_d: max_toporright + 5];
% Y = [min_leftorbottom - 5: h_d: max_toporright + 5];

X = [min_leftorbottom - 5: h_d: max_toporright + 5];
Y = [min_leftorbottom - 5: h_d: min_leftorbottom + 10 ];

U_i_PS = PS_incident_domain(xps1, yps1, X, Y, k) + PS_incident_domain(xps2, yps2, X, Y, k);

tic
[us1_0, x1, y1] = scattered_domain( G1, aj_1_r(:, 1), k, X, Y, N1);
figure(); pcolor(X, Y, real(us1_0)); shading interp; colorbar; title('scattered field no-rereflections ')
toc

UT_1_0 = U_i_PS - us1_0;

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

UT_2_1 = U_i_PS - us2_1;

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

for j = 2:3
    
    %solution on \Gamma_{1}
    [us1_r, ~, ~] = scattered_domain2( G1, G2, aj_1_r(:, j), aj_2_r(:, j-1), k, X, Y, N1, N2);
    UT_1_r = U_i_PS - us1_r;
    
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
    UT_2_r = U_i_PS - us2_r;
    
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


%% higher order scattered field
X_coordinates = [G1(1), G1(3), G2(1), G2(3)];
X_min = min(X_coordinates);
X_max = max(X_coordinates);
X_diff = abs(X_max - X_min);

Y_coordinates = [G1(2), G1(4), G2(2), G2(4)];
Y_min = min(Y_coordinates);
Y_max = max(Y_coordinates);
Y_diff = abs(Y_max - Y_min);

% maybe also change this part
N_diff = max(abs(Y_diff - X_diff)) + 2*pi;  % largest distance either x direction or y + 4*pi (2*pi in either direction)

N_d = 300;

min_leftorbottom = min(X_min, Y_min);
max_toporright = max(X_max, Y_max);
h_d = N_diff/N_d;
% need to change this bit so that it is closer zoomed in on the picture
% X = [min_leftorbottom - 5: h_d: max_toporright + 5];
% Y = [min_leftorbottom - 5: h_d: max_toporright + 5];

X = [min_leftorbottom - 5: h_d: max_toporright + 5];
Y = linspace(0, 15, length(X));
% Y = [0: h_d: max_toporright + 5 ];

U_i_PS = PS_incident_domain(xps1, yps1, X, Y, k) + PS_incident_domain(xps2, yps2, X, Y, k);

tic
[us2_highorder, x1, y1, x2, y2] = scattered_domain2( G1, G2, aj_1_r(:, end), aj_2_r(:, end), k, X, Y, N1, N2);
toc

U_T_highorder = U_i_PS - us2_highorder;
figure(); pcolor(X, Y, real(U_T_highorder)); shading interp; colormap(jet); colorbar; 
% title('total field 1 rereflections ')
hold on

Gamma_1 = plot(x1(length(x1)/2+1:end), y1(length(y1)/2+1:end));
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2(length(x2)/2+1:end), y2(length(y2)/2+1:end));
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];
% axis equal
% x_0 = [-dist: h_d: dist];
% y_0 = zeros(size(x_0));
% ground = plot(x_0, y_0);
% ground.LineWidth = 2;
% ground.Color = [0 0 0];

%% SPL computation
% first need to compute u_ref, the (complex valued) total field you would
% get if you were a distance R = 1m away from the source, with no barriers
% present.
u_ref = zeros(length(Y), length(X));
for xj = 1: length(X)
    
    for yj = 1:length(Y)
        
        u_ref( yj, xj) = 1i*besselh(0, k)/4;
        
        u_ref_norm( yj, xj) = norm(u_ref(yj, xj), 1);
        
        
        SPL(yj, xj) = 20*log10( norm(U_T_highorder(yj, xj), 1)/u_ref_norm(yj, xj) );

    end
    
end

% SPL = 20*log10( norm(U_T_highorder, 1)/u_ref_norm );
 
figure()
pcolor(X, Y, SPL); shading interp; colormap(jet); colorbar; 
hold on
Gamma_1 = plot(x1(length(x1)/2+1:end), y1(length(y1)/2+1:end));
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2(length(x2)/2+1:end), y2(length(y2)/2+1:end));
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];
% axis equal
