function [v_N_G2_r_1]  = multiple_scattering_2screen_step1(kwave, theta,...
    x2_col, y2_col, col_points2, vertices2, L2, x2nq2, y2nq2, nq2, tq2, ...
    h2, C1, C2, d, vertices1, L1, x1_t_outer, y1_t_outer, h1_outer, ...
    nq1_outer, x1_t_inner, y1_t_inner, h1_inner, nq1_inner, n2, v_N_G1_r0,...
    A2, VHNA2)
% In this function we will compute the coefficients for the 1st step. We
% will be evaluateing the integrals on the right hand side, there will be
% three components (f_{2}^{(1)}):
% incident wave
% -1i/4 \int_{0}^{L_{1}} H_{0}^{(1)}(k|x_{2}(s) - x_{1}(t)|)\phi_{1}^{0}(t)
% -1/4 \int_{0}^{L_{2}} H_{0}^{(1)}( k |s - t|) \Psi_{2}^{(1)}(t) dt
%
% \Psi_{1}^{0} = 2 \alpha_{1} dui/dn + ik/2 \int_{0}^{L_{1}} H_{1}^{1}(k | x_{2}(s) - x_{1}(t)|)
%                                        times (x_{2}(s) - x_{1}(t))/(
%                                        |x_{2}(s) - x_{1}(t)| ) \cdot n
%                                        \phi_{1}^{0} dt

% Inputs
%
% Outputs

% incident wave:
u_inc = incident(kwave, theta, x2_col, y2_col);

% Hankel function multiplied by phi_{1}^{(0)}
phi_1_0_int1 = v_N_G1_r0.eval(nq1_outer.', 1) + 2*duidn(vertices1, L1, kwave, d, nq1_outer.');
int_1 = midpoint_hankel_f_diff_screen(kwave, x2_col, y2_col, x1_t_outer, y1_t_outer, h1_outer, nq1_outer, phi_1_0_int1.' );

% Leading order behaviour
f_duidn = duidn(vertices2, L2, kwave, d, nq2);
LoB_1 = 2*PIM_int_hankel_f(kwave, col_points2, h2, nq2, f_duidn, tq2, C1, C2);

% now computing the double integral, doing slowly collocation point by
% collocation point

% inner integral: "collocation points" = integration nodes for \Gamma_{2}
%                   integration nodes = integration nodes are \Gamma_{1}
phi1_0_int_inner = v_N_G1_r0.eval(nq1_inner.', 1) + 2*duidn(vertices1, L1, kwave, d, nq1_inner.');

inner_int = midpoint_dphikdn_f_diff_screen(kwave, x2nq2, y2nq2, ...
    h1_inner, x1_t_inner, y1_t_inner, phi1_0_int_inner.', n2 );
% 
% inner_int = 1i*kwave*midpoint_dphikdn_f_diff_screen(kwave, x2nq2, y2nq2, ...
%     h1_inner, x1_t_inner, y1_t_inner, phi1_0_int_inner.', n2 )/2;

LoB_2 = PIM_int_hankel_f(kwave, col_points2, h2, nq2, inner_int.', tq2, C1, C2);

LoB = LoB_1 + LoB_2;
% putting altogether
RHS = u_inc - int_1 - LoB;

aj_2_1 = A2\RHS;

v_N_G2_r_1 = ProjectionFunction(aj_2_1, VHNA2);

