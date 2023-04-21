function [int, err, err_aver] = Psi_2_0_int_conv_test(kwave, d, Q_max, G1, vertices1, x2, y2, v_N_G1_r0, n2)
% In this function we will be testing the convergence of the integral
% \int_{0}^{L1}H_{1}^{1}(k|x2 - x1(t)|) (x2 -x1(t))/(|x2 - x1(t)|) \cdot
% n_{1}\phi1^0(t) dt
%
% For an increasing number of quadratuer points to see how it behaves, with
% a plot and errors I think
%
% Inputs:
%
% Outputs:
%
int = zeros(Q_max, length(x2));
for q = 1:Q_max

    Q = 2^q;  % number of quadrature points

    [x1, y1, ~, t1_mid, h1, ~, ~, L1] = ...
        discretisation_variables(G1, 1/Q, kwave);

    phi1_0 = v_N_G1_r0.eval(t1_mid.', 1) ...
        + 2*duidn(vertices1, L1, kwave, d, t1_mid.');

    int(q, :) = 1i*kwave*midpoint_dphikdn_f_diff_screen(kwave, x2, y2, ...
        h1, x1, y1, phi1_0.', n2)/2;

end

Q_true = 10*Q;

[x1, y1, ~, t1_mid, h1, ~, ~, L1] = ...
        discretisation_variables(G1, 1/Q_true, kwave);

 phi1_0 = v_N_G1_r0.eval(t1_mid.', 1) ...
        + 2*duidn(vertices1, L1, kwave, d, t1_mid.');

int_true = 1i*kwave*midpoint_dphikdn_f_diff_screen(kwave, x2, y2, ...
        h1, x1, y1, phi1_0.', n2)/2;
err = zeros(Q_max, length(x2));
err_aver = zeros(Q_max, 1);
for j = 1:Q_max
    err(j, :) = abs(int_true.' - int(q, :));

    err_aver(j, 1) = sum(err(j, :));
end
