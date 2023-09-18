% plotting for first few iterations in a general configuration

clear all
%add path for all code
addpath('./more_robust_attempt')
addpath('./Simon_code/')
% introducing the screens


% switch and things get interesting
G1 = [-2*pi, 2*pi, 2*pi, -2*pi];

G2 = [pi, 2*pi, 5*pi/2, pi]; 

C_wl= 1/20

k = 5;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% number of iterations for approximation and truth

R_it = 4;  % R_it = 13 (I think) for len =0.5. R_it > 50, len = 2;
R_true = 20;

[err_normG1_it, err_normG2_it, aj_1_r, aj_2_r, aj1_1step, aj2_1step] = it_conv_test_fixed_dof( G1, G2, k, C1, C2, theta, C_wl, R_it, R_true); 


% computing variables for G1:
[x1, y1, t1, h1, h1vector, N1, L1] = discretisation_variables(G1, C_wl, k);
% computing variables for G2:
[x2, y2, t2, h2, h2vector, N2, L2] = discretisation_variables(G2, C_wl, k);

%% Simons code set up
r1start = [G1(1), G1(2)];
r1end = [G1(3), G1(4)];

r2start = [G2(1), G2(2)];
r2end = [G2(3), G2(4)];


L1 = sqrt( (r1end(2) - r1start(2))^2 + ( r1end(1) - r1start(1) )^2 )
L2 = sqrt( (r2end(2) - r2start(2))^2 + ( r2end(1) - r2start(1) )^2 )

d = [0,-1]; % the direction of the incident wave
N1_SC = ceil(k*L1./(C_wl*2*pi)); N2_SC = ceil(k*L2./(C_wl*2*pi));% the number of boundary elements on screens 1 and 2

N_SC = max(N1_SC, N2_SC);
% N1_SC = N_SC
% N2_SC = N_SC
h1_SC = norm(r1end-r1start)/N1_SC; % the length of each element on screen 1
h2_SC = norm(r2end-r2start)/N2_SC; % the length of each element on screen 2
% Next calculate the x and y coordinates of the element midpoints on screen
% 1 and then on screen 2
x1_SC = r1start(1)+((1:N1_SC)-0.5)*(r1end(1)-r1start(1))/N1_SC;
y1_SC = r1start(2)+((1:N1_SC)-0.5)*(r1end(2)-r1start(2))/N1_SC;
x2_SC = r2start(1)+((1:N2_SC)-0.5)*(r2end(1)-r2start(1))/N2_SC;
y2_SC = r2start(2)+((1:N2_SC)-0.5)*(r2end(2)-r2start(2))/N2_SC;
h_SC = [h1_SC*ones(size(x1_SC)),h2_SC*ones(size(x2_SC))]; % the lengths of the elements
x_SC = [x1_SC,x2_SC]; % the x-coords of the element midpoints
y_SC = [y1_SC,y2_SC]; % and their y-coords
phi = bem2(x_SC,y_SC,h_SC,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
phi1 = phi(1:N1_SC); phi2 = phi(N1_SC+1:end); % Split phi into the solution on screens 1 and 2


%% computing error wrt r
for r = 1:R_true
    
    err_G1_L1norm(r, 1) = norm(phi1 - aj_1_r(:, r))/(norm(phi1));
    
    err_G2_L1norm(r, 1) = norm(phi2 - aj_2_r(:, r))/(norm(phi2));
    
end



%% plotting in domain
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


% solution in the domain

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

for j = 2:2
    
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

