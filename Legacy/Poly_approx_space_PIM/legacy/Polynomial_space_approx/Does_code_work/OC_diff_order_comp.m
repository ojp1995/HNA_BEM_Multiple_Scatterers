% this is a more explicit test to see if code differes depending on the way 
% the coordinates are input in the method 
%
% Oliver CODE
%
% We will run the code 4 times, once for the coordinates oriented as in
% other examples, once with the coordinates flipped for G1, another time
% with flipped for G2 and finally with both flipped.
% (By flipped I mean the start coordianted put in the end and the end put 
% in the start etc)
clear all

% control case
G1 = [0, 0, -2*pi, 2*pi];
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1

G2 = [4*pi, 0,  7*pi, 3*pi];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

%%% 
C_wl= 1/40

k = 5;  % wavenumber

N1 = ceil(k*L1./(C_wl*2*pi)) % number of itervals on G1
N2 = ceil(k*L2./(C_wl*2*pi)) % number of intervals on G2
 
theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi; 

Q = 1;

%%%% set up for other 3 cases:
G1A = [G1(3), G1(4), G1(1), G1(2)];
G2B = [G2(3), G2(4), G2(1), G2(2)];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  CONTROL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A1_con, x_col1 A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
[A2_con, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);
[S21_con, col_points_1] = S21_op(G1, G2, k, N1, N2);
[S12_con, col_points_2] = S21_op(G2, G1, k, N2, N1);

u_i_1_0_con = PW_incident(k, theta, G1, x_col1);
aj_1_0_con = A1_con\u_i_1_0_con.';

 % step 1
u_i_2_1_con =  PW_incident(k, theta, G2, x_col2).' - S21_con*coeff_2soln_midpoint(aj_1_0_con, L1, N1-1, N1);
aj_2_1_con = A2_con\u_i_2_1_con;

% step 2
u_i_1_2_con =  PW_incident(k, theta, G1, x_col1).' - S12_con*coeff_2soln_midpoint(aj_2_1_con, L2, N2-1, N2);
aj_1_2_con = A1_con\u_i_1_2_con;

% step 3
u_i_2_3_con =  PW_incident(k, theta, G2, x_col2).' - S21_con*coeff_2soln_midpoint(aj_1_2_con, L1, N1-1, N1);
aj_2_3_con = A2_con\u_i_2_3_con;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CASE A: G1 coordinates swapped %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A1_A, x_col1 A1_singA, A1_smoothA] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
[A2_A, x_col2, A2_singA, A2_smoothA] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);
[S21_A, col_points_1] = S21_op(G1A, G2, k, N1, N2);
[S12_A, col_points_2] = S21_op(G2, G1A, k, N2, N1);

u_i_1_0_A = PW_incident(k, theta, G1A, x_col1);
aj_1_0_A = A1_A\u_i_1_0_A.';

 % step 1
u_i_2_1_A =  PW_incident(k, theta, G2, x_col2).' - fliplr(S21_A)*coeff_2soln_midpoint(flipud(aj_1_0_A), L1, N1-1, N1);
aj_2_1_A = A2_A\u_i_2_1_A;

% step 2
u_i_1_2_A =  PW_incident(k, theta, G1A, x_col1).' - S12_A*coeff_2soln_midpoint(aj_2_1_A, L2, N2-1, N2);
aj_1_2_A = A1_A\u_i_1_2_A;

% step 3
u_i_2_3_A =  PW_incident(k, theta, G2, x_col2).' - fliplr(S21_A)*coeff_2soln_midpoint(flipud(aj_1_2_A), L1, N1-1, N1);
aj_2_3_A = A2_A\u_i_2_3_A;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Comparison to control case%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(x_col1/L1, real(aj_1_0_con - aj_1_0_A), 'DisplayName', 'difference between case A and control')
hold on
plot(x_col1/L1, real(aj_1_0_con), 'DisplayName', 'control')
plot(x_col1/L1, real(aj_1_0_A), 'DisplayName', 'Case A')
plot(x_col1/L1, real(aj_1_0_con - flipud(aj_1_0_A)), 'DisplayName', 'difference between case A and control')
legend show
xlabel('x/L')
ylabel('$\phi_{1}^{0}$ ')
title('Approximations of $\phi_{1}^{(0)}$, CASE A')
%%
figure()
plot(x_col2/L2, real(aj_2_1_con - aj_2_1_A), 'DisplayName', 'difference between case A and control')
hold on
% plot(x_col2/L2, real(aj_2_1_con), 'DisplayName', 'control')
% plot(x_col2/L2, real(aj_2_1_A), 'DisplayName', 'Case A')
legend show
xlabel('x/L')
ylabel('$\phi_{2}^{1}$ ')
title('Approximations of $\phi_{2}^{(1)}$, CASE A')

% this is massib=vely different, so what is different between the two
% cases?
% 1. S21 matrix
% 2. Plane wave being incident on G2? --> It cannot be, it is the same code

S21_A_vecdiff = S21_con*coeff_2soln_midpoint(aj_1_0_con, L1, N1-1, N1) - S21_A*coeff_2soln_midpoint(aj_1_0_A, L1, N1-1, N1);
S21_A_vecdiff_A = S21_con*aj_1_0_con - S21_A*aj_1_0_A;
S21_A_matdiff = S21_con - S21_A; 
figure()
surf(real(S21_A_matdiff))
% is it a rotation?
% figure()
% surf(real(S21_con - flipud(S21_A)))
% title('flipud down diff')

figure()
surf(real(S21_con - fliplr(S21_A)))
title('fliplr diff')

% We need to flip the matrix lr.
 % step 1
u_i_2_1_A_lr =  PW_incident(k, theta, G2, x_col2).' - fliplr(S21_A)*flipud(aj_1_0_A);
aj_2_1_A_lr = A2_A\u_i_2_1_A_lr;

figure()
plot(x_col2/L2, real(aj_2_1_con - aj_2_1_A_lr), 'DisplayName', 'difference between case A and control')
hold on
title('flipped S21 lr, aj_1_0 also flipud')

%%
figure()
plot(x_col1/L1, real(aj_1_2_con - flipud(aj_1_2_A)), 'DisplayName', 'flipud difference between case A and control')
hold on
plot(x_col1/L1, real(aj_1_2_con), 'DisplayName', 'control')
plot(x_col1/L1, real(aj_1_2_A), 'DisplayName', 'Case A')
plot(x_col1/L1, real(aj_1_2_con - aj_1_2_A), 'DisplayName', 'difference between case A and control')

legend show
xlabel('x/L')
ylabel('$\phi_{1}^{0}$ ')
title('Approximations of $\phi_{1}^{(2)}$, CASE A')

figure()
plot(x_col2/L2, real(aj_2_3_con - aj_2_3_A), 'DisplayName', 'difference between case A and control')
hold on
plot(x_col2/L2, real(aj_2_3_con), 'DisplayName', 'control')
plot(x_col2/L2, real(aj_2_3_A), 'DisplayName', 'Case A')
legend show
xlabel('x/L')
ylabel('$\phi_{2}^{3}$ ')
title('Approximations of $\phi_{2}^{(3)}$, CASE A')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  CASE B: G2 coordinates swapped %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A1_B, x_col1 A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
[A2_B, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);
[S21_B, col_points_1] = S21_op(G1, G2B, k, N1, N2);
[S12_B, col_points_2] = S21_op(G2B, G1, k, N2, N1);

%%
%%%% Comparison of the S21 operators
figure(); surf(real(S21_B)); shading interp; title('S21B operator')
figure(); surf(real(S21_con)); shading interp; title('S21 control operator')
figure(); surf(real(S21_con - S21_B)); shading interp; title('S21 control - S21B')
figure(); surf(real(S21_con - flipud(S21_B))); shading interp; title('S21 control - flipud(S21B)')
figure(); surf(real(S21_con - fliplr(S21_B))); shading interp; title('S21 control - fliplr(S21B)')

%%%% comparison of A2 matrix
figure(); surf(real(A2_B)); shading interp; title('A2B matrix')
figure(); surf(real(A2_con)); shading interp; title('A2 control matrix')
figure(); surf(real(A2_con - A2_B)); shading interp; title('Difference, A2con - A2')
%%

%%%%%% Doesn't work!
% % % % % % % % flipping S21_B operator
% S21_B = flipud(S21_B);


u_i_1_0_B = PW_incident(k, theta, G1, x_col1);
aj_1_0_B = A1_B\u_i_1_0_B.';

 % step 1
u_i_2_1_B =  PW_incident(k, theta, G2B, x_col2).' - S21_B*coeff_2soln_midpoint(aj_1_0_B, L1, N1-1, N1);
% u_i_2_1_B =  flipud(PW_incident(k, theta, G2B, x_col2).') - S21_B*coeff_2soln_midpoint(aj_1_0_B, L1, N1-1, N1);
aj_2_1_B = A2_B\u_i_2_1_B;

% step 2
u_i_1_2_B =  PW_incident(k, theta, G1, x_col1).' - S12_B*coeff_2soln_midpoint(aj_2_1_B, L2, N2-1, N2);
aj_1_2_B = A1_B\u_i_1_2_B;

% step 3
u_i_2_3_B =  PW_incident(k, theta, G2B, x_col2).' - S21_B*coeff_2soln_midpoint(aj_1_2_B, L1, N1-1, N1);
aj_2_3_B = A2_B\u_i_2_3_B;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Comparison to control case%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conmparison between RHS terms
figure(); plot(x_col2/L2, real(PW_incident(k, theta, G2B, x_col2) - PW_incident(k, theta, G2, x_col2))); title('Comparison of plane waves')
figure(); plot(x_col2/L2, real(PW_incident(k, theta, G2B, x_col2) - fliplr(PW_incident(k, theta, G2, x_col2)))); title('Comparison of plane waves, fliplr')

% comparisons to coefficients
figure()
plot(x_col1/L1, real(aj_1_0_con - aj_1_0_B), 'DisplayName', 'difference between case B and control')
hold on
plot(x_col1/L1, real(aj_1_0_con), 'DisplayName', 'control')
plot(x_col1/L1, real(aj_1_0_B), 'DisplayName', 'Case B')
legend show
xlabel('x/L')
ylabel('$\phi_{1}^{0}$ ')
title('Approximations of $\phi_{1}^{(0)}$, CASE B')

figure()
plot(x_col2/L2, real(aj_2_1_con - aj_2_1_B), 'DisplayName', 'difference between case B and control')
hold on
plot(x_col2/L2, real(aj_2_1_con), 'DisplayName', 'control')
plot(x_col2/L2, real(aj_2_1_B), 'DisplayName', 'Case B')
plot(x_col2/L2, real(flipud(aj_2_1_B)), 'DisplayName', 'Flipped Case B', 'LineStyle', '-.')
plot(x_col2/L2, real(aj_2_1_con - flipud(aj_2_1_B)), 'DisplayName', 'Diff, Case B flipped', 'LineStyle', ':')
legend show
xlabel('x/L')
ylabel('$\phi_{2}^{1}$ ')
title('Approximations of $\phi_{2}^{(1)}$, CASE B')

figure()
plot(x_col1/L1, real(aj_1_2_con - aj_1_2_B), 'DisplayName', 'difference between case B and control')
hold on
plot(x_col1/L1, real(aj_1_2_con), 'DisplayName', 'control')
plot(x_col1/L1, real(aj_1_2_B), 'DisplayName', 'Case B')
legend show
xlabel('x/L')
ylabel('$\phi_{1}^{2}$ ')
title('Approximations of $\phi_{1}^{(2)}$, CASE B')

figure()
plot(x_col2/L2, real(aj_2_3_con - aj_2_3_B), 'DisplayName', 'difference between case B and control')
hold on
plot(x_col2/L2, real(aj_2_3_con), 'DisplayName', 'control')
plot(x_col2/L2, real(aj_2_3_B), 'DisplayName', 'Case B')
legend show
xlabel('x/L')
ylabel('$\phi_{2}^{3}$ ')
title('Approximations of $\phi_{2}^{(3)}$, CASE B')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CASE C: G1 and G2 coordinates swapped %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A1_C, x_col1 A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
[A2_C, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);
[S21_C, col_points_1] = S21_op(G1A, G2B, k, N1, N2);
[S12_C, col_points_2] = S21_op(G2B, G1A, k, N2, N1);

u_i_1_0_C = PW_incident(k, theta, G1A, x_col1);
aj_1_0_C = A1_C\u_i_1_0_C.';

 % step 1
u_i_2_1_C =  PW_incident(k, theta, G2B, x_col2).' - S21_C*coeff_2soln_midpoint(aj_1_0_C, L1, N1-1, N1);
aj_2_1_C = A2_C\u_i_2_1_C;

% step 2
u_i_1_2_C =  PW_incident(k, theta, G1A, x_col1).' - S12_C*coeff_2soln_midpoint(aj_2_1_C, L2, N2-1, N2);
aj_1_2_C = A1_C\u_i_1_2_C;

% step 3
u_i_2_3_C =  PW_incident(k, theta, G2B, x_col2).' - S21_C*coeff_2soln_midpoint(aj_1_2_C, L1, N1-1, N1);
aj_2_3_C = A2_C\u_i_2_3_C;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Comparison to control case%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(x_col1/L1, real(aj_1_0_con - aj_1_0_C), 'DisplayName', 'difference between case C and control')
hold on
plot(x_col1/L1, real(aj_1_0_con), 'DisplayName', 'control')
plot(x_col1/L1, real(aj_1_0_C), 'DisplayName', 'Case C')
legend show
xlabel('x/L')
ylabel('$\phi_{1}^{0}$ ')
title('Approximations of $\phi_{1}^{(0)}$, CASE C')

figure()
plot(x_col2/L2, real(aj_2_1_con - aj_2_1_C), 'DisplayName', 'difference between case C and control')
hold on
plot(x_col2/L2, real(aj_2_1_con), 'DisplayName', 'control')
plot(x_col2/L2, real(aj_2_1_C), 'DisplayName', 'Case C')
legend show
xlabel('x/L')
ylabel('$\phi_{2}^{1}$ ')
title('Approximations of $\phi_{2}^{(1)}$, CASE C')

figure()
plot(x_col1/L1, real(aj_1_2_con - aj_1_2_C), 'DisplayName', 'difference between case C and control')
hold on
plot(x_col1/L1, real(aj_1_2_con), 'DisplayName', 'control')
plot(x_col1/L1, real(aj_1_2_C), 'DisplayName', 'Case C')
legend show
xlabel('x/L')
ylabel('$\phi_{1}^{2}$ ')
title('Approximations of $\phi_{1}^{(2)}$, CASE C')

figure()
plot(x_col2/L2, real(aj_2_3_con - aj_2_3_C), 'DisplayName', 'difference between case C and control')
hold on
plot(x_col2/L2, real(aj_2_3_con), 'DisplayName', 'control')
plot(x_col2/L2, real(aj_2_3_C), 'DisplayName', 'Case C')
legend show
xlabel('x/L')
ylabel('$\phi_{2}^{3}$ ')
title('Approximations of $\phi_{2}^{(3)}$, CASE C')



