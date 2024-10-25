function [v_N_G1_r0, RHS] = multiple_Scattering_2screen_step0(kwave, theta, ...
    d, vertices1, L1, col_points1, x1_col, y1_col, h, nq, tq, A, VHNA1,...
    C1, C2, alpha1)
% In this function we will be performing the computations for step 0.
% ...
% Input:
% kwave - wavenumber
% theta - angle of incident wave
% d - direction of incident wave
% vertices1 - matrix of coordinates of screen
% L1 - length of \Gamma_{1}
% col_points1 - parameterised collocation points
% (x1_col, y1_col) - cartesian collocation points
% nq are the quadrature nodes
% tq is the discretisation mesh for the quadrature.
% A - Left hand side Matrix
% VHNA1 - object we are copying structure of/borrowing
% alpha1
%
% Output:
% v_N_G1_r0 - object which contains the coefficients
% RHS - RHS vec just incase


% Need to compute the right hand side, two steps, leading order term
% integrated and then also the incident wave

u_inc = incident(kwave, theta, x1_col, y1_col);

f_duidn = duidn(vertices1, L1, kwave, d, nq);
LoB = 2*alpha1*PIM_int_hankel_f(kwave, col_points1, h, nq, f_duidn, tq, C1, C2);

RHS = u_inc - LoB;

aj = A\RHS;

v_N_G1_r0 = ProjectionFunction(aj, VHNA1);



